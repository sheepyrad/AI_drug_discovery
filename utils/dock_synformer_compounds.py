#!/usr/bin/env python
"""
Module for docking top compounds from Synformer retrosynthesis results.

Provides a function to dock the top synthesizable compounds from Synformer
retrosynthesis analysis, using the same parameters as the main pipeline.
"""

import logging
from pathlib import Path
import pandas as pd
import shutil
import os

# Import the redocking function
from .redocking import redock_compound, vfu_dir

logger = logging.getLogger(__name__)

def dock_synformer_compounds(retro_dir, round_no, receptor, out_dir, top_n=5, 
                             program_choice="qvina", scoring_function="nnscore2",
                             center=(114.817, 75.602, 82.416), box_size=(38, 70, 58),
                             exhaustiveness=10, is_selfies=False, is_peptide=False,
                             log_callback=None):
    """
    Dock the top synthesizable compounds from Synformer results.
    
    Args:
        retro_dir: Path to the retrosynthesis results directory
        round_no: Round number (0-based)
        receptor: Filename of the receptor file (should be in ./input directory)
        out_dir: Output directory for saving results
        top_n: Number of top compounds to dock (default: 5)
        program_choice: Docking program choice (default: "qvina")
        scoring_function: Scoring function (default: "nnscore2")
        center: X Y Z coordinates of the docking center
        box_size: Box sizes in X Y Z
        exhaustiveness: Exhaustiveness for docking
        is_selfies: Use SELFIES in redocking (default: False)
        is_peptide: Flag if the ligand is a peptide (default: False)
        log_callback: Function to call for logging (default: None)
    
    Returns:
        Path to the directory with docking results
    """
    # Use the provided logging callback or the default logger
    log = log_callback or logger.info
    
    # Create full paths
    retro_dir = Path(retro_dir)
    out_dir = Path(out_dir)
    # No need to convert receptor to absolute path - it should be a filename
    
    # Ensure output directory exists
    out_dir.mkdir(exist_ok=True, parents=True)
    
    log(f"Docking top {top_n} compounds from Synformer results for round {round_no+1}")
    
    # Create directory for Synformer docking results
    synformer_dock_dir = out_dir / f"synformer_docking_round_{round_no+1}"
    synformer_dock_dir.mkdir(exist_ok=True)
    
    # Path to round-specific retrosynthesis results
    round_retro_dir = retro_dir / f"round_{round_no+1}"
    
    if not round_retro_dir.exists():
        logger.error(f"Retrosynthesis directory not found: {round_retro_dir}")
        return None
    
    # Get all CSV files with retrosynthesis results (except the summary file)
    # Improved filtering to avoid picking up summary files
    retro_files = [f for f in round_retro_dir.glob("top*.csv") if not f.name.endswith("_summary.csv")]
    
    if not retro_files:
        log(f"No retrosynthesis results found in {round_retro_dir}")
        return synformer_dock_dir
    
    # Prepare docking parameters
    center_x, center_y, center_z = center
    size_x, size_y, size_z = box_size
    
    redock_params = (
        program_choice,
        scoring_function,
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z,
        exhaustiveness,
        is_selfies,
        is_peptide
    )
    
    results = []
    
    # Process each file (already sorted by rank because of filename pattern top1_, top2_, etc.)
    for file_idx, retro_file in enumerate(sorted(retro_files)[:top_n]):
        log(f"Processing retrosynthesis file: {retro_file.name}")
        
        try:
            # Read the CSV file
            df = pd.read_csv(retro_file)
            
            if df.empty:
                log(f"Empty retrosynthesis file: {retro_file}")
                continue
            
            # Process top 5 synthetic products from each file (or fewer if the file has fewer rows)
            product_count = min(5, len(df))
            
            for row_idx in range(product_count):
                # KEY CHANGE: Prioritize the 'smiles' column over 'target'
                # The 'smiles' column contains proposed synthetic products
                if 'smiles' in df.columns:
                    smiles = df['smiles'].iloc[row_idx]
                else:
                    # If column structure is different, take the first SMILES-looking column
                    for col in df.columns:
                        if any(c for c in df[col].iloc[row_idx] if c in '()=#@'):
                            smiles = df[col].iloc[row_idx]
                            break
                    else:
                        log(f"Could not find SMILES in {retro_file}, row {row_idx}")
                        continue
                
                # If the product has a score, use it for identification
                score_text = ""
                if 'score' in df.columns:
                    score = df['score'].iloc[row_idx]
                    score_text = f"_score{score:.3f}"
                
                # Generate a unique ID for this synthetic product
                orig_compound = retro_file.stem.replace("top", "")
                cid = f"synth{file_idx+1}_{row_idx+1}{score_text}_{orig_compound}"
                
                log(f"Docking synthetic product {row_idx+1} from {retro_file.name}: {smiles}")
                
                # Perform redocking
                pose_out, rescored = redock_compound(
                    cid,
                    smiles,
                    redock_params,
                    receptor=receptor,
                    log_callback=log
                )
                
                # Store results
                results.append({
                    "compound_id": cid,
                    "smiles": smiles,
                    "pose_pred_out": pose_out,
                    "re_scored_values": rescored
                })
                
                # Copy poses for this compound
                vfu_outputs_dir = Path(vfu_dir) / "outputs"
                compound_poses_dir = synformer_dock_dir / f"synth_compound_{cid}"
                compound_poses_dir.mkdir(exist_ok=True)
                
                # Copy all files from VFU outputs to the compound-specific directory
                for file_path in vfu_outputs_dir.glob("*"):
                    if file_path.is_file():
                        shutil.copy2(file_path, compound_poses_dir)
                    elif file_path.is_dir():
                        dest_dir = compound_poses_dir / file_path.name
                        if dest_dir.exists():
                            shutil.rmtree(dest_dir)
                        shutil.copytree(file_path, dest_dir)
                
                log(f"Copied poses for compound {cid} to {compound_poses_dir}")
            
        except Exception as e:
            logger.error(f"Error processing {retro_file}: {e}")
    
    # Save results to CSV
    if results:
        results_df = pd.DataFrame(results)
        results_csv = synformer_dock_dir / f"synformer_docking_results_round_{round_no+1}.csv"
        results_df.to_csv(results_csv, index=False)
        log(f"Synformer docking results saved to: {results_csv}")
    
    return synformer_dock_dir 