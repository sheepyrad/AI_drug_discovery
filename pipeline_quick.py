#!/usr/bin/env python
"""
Quick pipeline for de novo drug discovery.

Streamlined workflow:
  1. Ligand Generation (same as previous pipeline)
  2. Convert SDF to SMILES strings
  3. Run retrosynthesis (Synformer) on each compound
  4. Extract top 5 variants from retrosynthesis results
  5. Apply medchem filtering to the variants
  6. Redock filtered variants to the receptor

The quick pipeline skips energy minimization, posebuster evaluation, 
and other intermediate steps to provide faster results.
"""

import os
import sys
import time
from pathlib import Path
import logging
import pandas as pd
import shutil
from rdkit import Chem
import tempfile
import multiprocessing
from multiprocessing import Process, Queue
from contextlib import contextmanager
import queue
import threading

# Set up thread-safe logging with QueueHandler and QueueListener
class ThreadSafeRotatingFileHandler(logging.FileHandler):
    """A file handler that is thread-safe and properly handles file descriptor issues"""
    def __init__(self, filename, mode='a', encoding=None, delay=False):
        super().__init__(filename, mode, encoding, delay)
        self._lock = threading.RLock()
        
    def emit(self, record):
        with self._lock:
            try:
                super().emit(record)
            except Exception:
                # If there's an error with the file descriptor, reopen the file
                try:
                    self.close()
                    self.stream = self._open()
                    super().emit(record)
                except Exception:
                    self.handleError(record)

# Create a queue for thread-safe logging
log_queue = queue.Queue(-1)  # No limit on size

# Configure the root logger with queue handler
root_logger = logging.getLogger()
root_logger.setLevel(logging.INFO)

# Create handlers
console_handler = logging.StreamHandler(sys.stdout)
file_handler = ThreadSafeRotatingFileHandler("quick_pipeline.log", mode="w")

# Set formatter
formatter = logging.Formatter("[%(asctime)s] %(name)s - %(levelname)s - %(message)s")
console_handler.setFormatter(formatter)
file_handler.setFormatter(formatter)

# Add handlers to root logger
root_logger.addHandler(console_handler)
root_logger.addHandler(file_handler)

# Get logger for this module
logger = logging.getLogger(__name__)
logger.info("Starting quick pipeline...")

# Import functions from the submodules
from ligand_generation import run_ligand_generation
from redocking import redock_compound, vfu_dir
from retrosynformer import run_retrosynthesis
from medchem_filter import generative_filter, filter_compounds

def extract_smiles_from_sdf(sdf_file):
    """
    Extract SMILES strings from an SDF file.
    
    Args:
        sdf_file: Path to the SDF file
        
    Returns:
        List of dicts with compound_id and smiles
    """
    compounds = []
    try:
        # Read molecules from SDF
        supplier = Chem.SDMolSupplier(str(sdf_file))
        
        for idx, mol in enumerate(supplier):
            if mol is not None:
                # Generate compound ID based on the file name and index
                compound_id = f"{sdf_file.stem}_mol_{idx+1}"
                smiles = Chem.MolToSmiles(mol)
                compounds.append({"compound_id": compound_id, "smiles": smiles})
    except Exception as e:
        logger.error(f"Error extracting SMILES from SDF: {e}")
    
    logger.info(f"Extracted {len(compounds)} SMILES strings from {sdf_file}")
    return compounds

def extract_variants_from_retrosynthesis(retro_result_file, max_variants=5):
    """
    Extract synthetic variants from a retrosynthesis result file.
    
    Args:
        retro_result_file: Path to the CSV file with retrosynthesis results
        max_variants: Maximum number of variants to extract
        
    Returns:
        List of dicts with variant_id and smiles
    """
    variants = []
    try:
        # Read retrosynthesis results
        df = pd.read_csv(retro_result_file)
        
        if df.empty:
            logger.warning(f"Empty retrosynthesis results file: {retro_result_file}")
            return variants
        
        # Find the SMILES column (should be 'smiles' but check alternatives)
        smiles_col = None
        if 'smiles' in df.columns:
            smiles_col = 'smiles'
        else:
            # Look for a column that might contain SMILES strings
            for col in df.columns:
                if df[col].dtype == 'object' and any(c for c in df[col].iloc[0] if c in '()=#@'):
                    smiles_col = col
                    break
        
        if smiles_col is None:
            logger.error(f"No SMILES column found in {retro_result_file}")
            return variants
        
        # Try to sort by a score column if it exists to get the best variants
        if 'score' in df.columns:
            # Determine if higher or lower score is better (assuming higher is better by default)
            # NOTE: If your scores are like energy where lower is better, reverse the sort
            df = df.sort_values('score', ascending=False)
            logger.info(f"Sorting variants by 'score' column (higher is better)")
        
        # Extract up to max_variants
        for idx, row in df.head(max_variants).iterrows():
            smiles = row[smiles_col]
            # Get score if available for labeling
            score_text = ""
            if 'score' in df.columns:
                score = row['score']
                if not pd.isna(score):
                    score_text = f"_score{score:.3f}"
            
            # Create variant ID based on the original compound and variant number
            parent_id = retro_result_file.stem
            variant_id = f"{parent_id}_variant_{idx+1}{score_text}"
            
            variants.append({
                "variant_id": variant_id, 
                "smiles": smiles, 
                "parent_id": parent_id,
                "score": float(row['score']) if 'score' in df.columns and not pd.isna(row['score']) else None
            })
    
    except Exception as e:
        logger.error(f"Error extracting variants from {retro_result_file}: {e}")
    
    logger.info(f"Extracted {len(variants)} variants from {retro_result_file}")
    return variants

def smiles_to_sdf(smiles_list, output_file):
    """
    Convert a list of SMILES strings to an SDF file.
    
    Args:
        smiles_list: List of dicts with 'smiles' and ID ('variant_id' or other id field)
        output_file: Path to save the SDF file
        
    Returns:
        Path to the created SDF file or None if failed
    """
    try:
        writer = Chem.SDWriter(str(output_file))
        
        for item in smiles_list:
            smiles = item['smiles']
            
            # Get the ID (could be compound_id, variant_id, etc.)
            id_field = next((k for k in item.keys() if k.endswith('_id')), None)
            id_value = item.get(id_field, "unknown")
            
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Add properties to the molecule
                mol.SetProp("_Name", id_value)
                mol.SetProp("SMILES", smiles)
                
                # Add barcode and tracking information
                if 'barcode' in item:
                    mol.SetProp("BARCODE", item['barcode'])
                if 'parent_id' in item:
                    mol.SetProp("PARENT_ID", item['parent_id'])
                if 'source_compound' in item:
                    mol.SetProp("SOURCE_COMPOUND", item['source_compound'])
                if 'generation' in item:
                    mol.SetProp("GENERATION", item['generation'])
                
                writer.write(mol)
            else:
                logger.warning(f"Could not convert SMILES to molecule: {smiles}")
        
        writer.close()
        logger.info(f"Created SDF file with {len(smiles_list)} molecules: {output_file}")
        return output_file
    
    except Exception as e:
        logger.error(f"Error converting SMILES to SDF: {e}")
        return None

def apply_medchem_filtering_to_variants(variants, output_dir):
    """
    Apply medchem filtering to variant SMILES.
    
    Args:
        variants: List of variant dictionaries with SMILES
        output_dir: Directory to save temporary files and results
        
    Returns:
        List of filtered variants
    """
    # Create a temporary SDF file with all variants
    temp_sdf = output_dir / "temp_variants_for_filtering.sdf"
    smiles_to_sdf(variants, temp_sdf)
    
    try:
        # Apply MedChem filtering
        logger.info(f"Applying MedChem filtering to {len(variants)} variants...")
        filtered_df = generative_filter(temp_sdf, output_folder=output_dir)
        
        if filtered_df.empty:
            logger.warning("No variants passed MedChem filtering")
            return []
        
        # Create a list of filtered variants
        filtered_variants = []
        
        # Log the column names to help with debugging
        logger.info(f"Filtered dataframe columns: {filtered_df.columns.tolist()}")
        
        # Determine which column contains the variant ID
        id_column = None
        if 'Name' in filtered_df.columns:
            id_column = 'Name'
        elif 'compound_id' in filtered_df.columns:
            id_column = 'compound_id'
        else:
            # If neither column exists, log the issue and return empty list
            logger.error(f"Could not find Name or compound_id column in filtered results. Available columns: {filtered_df.columns.tolist()}")
            return []
            
        logger.info(f"Using '{id_column}' column to match filtered variants")
        
        # Create a dictionary for faster lookup of original variants
        variant_lookup = {
            v.get('variant_id', ''): v for v in variants
        }
        
        # Also try matching by barcode if exists
        barcode_lookup = {}
        for v in variants:
            if 'barcode' in v:
                barcode_lookup[v['barcode']] = v
        
        for _, row in filtered_df.iterrows():
            # Find the original variant info using the identified column
            variant_id = row[id_column]
            
            # Log the variant ID for debugging
            logger.info(f"Looking for variant with ID: {variant_id}")
            
            # First, try direct lookup by variant_id
            original_variant = variant_lookup.get(variant_id)
            
            # If not found, try other methods
            if not original_variant:
                # Try matching by barcode if available
                if 'BARCODE' in filtered_df.columns:
                    barcode = row['BARCODE']
                    original_variant = barcode_lookup.get(barcode)
                
                # If still not found, try iterating through variants
                if not original_variant:
                    for v in variants:
                        # Check if any ID field matches
                        for key, value in v.items():
                            if key.endswith('_id') and value == variant_id:
                                original_variant = v
                                break
                        if original_variant:
                            break
            
            if original_variant:
                logger.info(f"Found matching variant: {original_variant.get('variant_id', 'unknown')} (Barcode: {original_variant.get('barcode', 'unknown')})")
                filtered_variants.append(original_variant)
            else:
                logger.warning(f"Could not find original variant for {variant_id}")
        
        logger.info(f"After MedChem filtering, {len(filtered_variants)} variants remain")
        return filtered_variants
    
    except Exception as e:
        logger.error(f"Error during MedChem filtering: {e}")
        return []
    
    finally:
        # Clean up temporary files
        if temp_sdf.exists():
            temp_sdf.unlink()

def generate_tracking_report(compounds, variants, redock_results, out_dir):
    """
    Generate a comprehensive report showing the lineage of all compounds.
    
    Args:
        compounds: List of original generated compounds
        variants: List of variants from retrosynthesis
        redock_results: List of final docking results
        out_dir: Directory to save the report
    """
    try:
        # Create a comprehensive dataframe tracing all compounds
        report_rows = []
        
        # Add original compounds
        for comp in compounds:
            row = {
                'compound_id': comp.get('compound_id', 'unknown'),
                'barcode': comp.get('barcode', 'unknown'),
                'generation': comp.get('generation', '1'),
                'smiles': comp.get('smiles', ''),
                'parent_id': 'NONE',
                'status': 'GENERATED',
                'source': 'AI_GENERATION'
            }
            report_rows.append(row)
        
        # Add all variants
        for var in variants:
            row = {
                'compound_id': var.get('variant_id', 'unknown'),
                'barcode': var.get('barcode', 'unknown'),
                'generation': var.get('generation', '2'),
                'smiles': var.get('smiles', ''),
                'parent_id': var.get('source_compound', 'unknown'),
                'parent_barcode': next((c.get('barcode') for c in compounds if c.get('compound_id') == var.get('source_compound')), 'unknown'),
                'status': 'SYNTHETIZED',
                'source': 'RETROSYNTHESIS',
                'score': var.get('score', None)
            }
            report_rows.append(row)
        
        # Update status for successfully docked variants
        docked_barcodes = {result.get('barcode') for result in redock_results}
        for row in report_rows:
            if row['barcode'] in docked_barcodes:
                row['status'] = 'DOCKED'
                
                # Find the docking score
                for result in redock_results:
                    if result.get('barcode') == row['barcode']:
                        # Extract best docking score and pose from pose_pred_out
                        if 'pose_pred_out' in result and result['pose_pred_out'] is not None:
                            try:
                                import ast
                                # Check if the pose_pred_out is already a dictionary
                                if isinstance(result['pose_pred_out'], dict):
                                    pose_dict = result['pose_pred_out']
                                else:
                                    # Try to safely evaluate the string representation
                                    pose_dict = ast.literal_eval(result['pose_pred_out'])
                                
                                # Find the pose with the best (lowest) score
                                best_score = float('inf')
                                best_pose = None
                                
                                for pose_name, pose_data in pose_dict.items():
                                    # The pose_data is a list where the first element is a list of scores
                                    # and the second element is the path to the pose file
                                    if pose_data and isinstance(pose_data[0], list) and pose_data[0]:
                                        score = pose_data[0][0]  # Get the first score from the first list
                                        if score < best_score:
                                            best_score = score
                                            best_pose = pose_name
                                
                                if best_pose:
                                    row['docking_score'] = best_score
                                    row['best_pose'] = best_pose.replace('.pdbqt', '')
                            except Exception as e:
                                logger.warning(f"Error extracting pose data for {row['barcode']}: {e}")
                                # Try to extract the best score directly from the string if parsing failed
                                try:
                                    if isinstance(result['pose_pred_out'], str):
                                        # Look for patterns like [[-6.2, ...]] to extract the best score
                                        import re
                                        score_match = re.search(r'\[\[([-\d\.]+)', result['pose_pred_out'])
                                        if score_match:
                                            row['docking_score'] = float(score_match.group(1))
                                            logger.info(f"Extracted fallback score {row['docking_score']} for {row['barcode']}")
                                except Exception as fallback_error:
                                    logger.warning(f"Fallback score extraction failed for {row['barcode']}: {fallback_error}")
                        # Fallback to re_scored_values if pose_pred_out processing fails
                        elif 're_scored_values' in result and result['re_scored_values'] is not None:
                            try:
                                docking_score = result['re_scored_values'].split(',')[0]
                                row['docking_score'] = float(docking_score)
                            except Exception as e:
                                logger.warning(f"Error extracting re_scored_values for {row['barcode']}: {e}")
        
        # Create the report DataFrame
        report_df = pd.DataFrame(report_rows)
        
        # Save to CSV
        report_file = out_dir / "compound_tracking_report.csv"
        report_df.to_csv(report_file, index=False)
        logger.info(f"Compound tracking report saved to: {report_file}")
        
        # Return the dataframe for further use if needed
        return report_df
    
    except Exception as e:
        logger.error(f"Error generating tracking report: {e}")
        return None

# Helper function to run retrosynthesis with timeout
def run_retrosynthesis_with_timeout(smiles, output_path, timeout=300):
    """
    Run retrosynthesis with a timeout using multiprocessing.
    
    Args:
        smiles: SMILES string of the compound
        output_path: Path to save the output CSV
        timeout: Maximum time in seconds to wait (default: 300 seconds / 5 minutes)
        
    Returns:
        True if successful, False if failed or timed out
    """
    def worker(smiles, output_path, queue):
        result = run_retrosynthesis(smiles, output_path)
        queue.put(result)
    
    # Create a queue to get the result
    queue = Queue()
    
    # Create and start the process
    process = Process(target=worker, args=(smiles, output_path, queue))
    process.start()
    
    # Wait for the specified timeout
    process.join(timeout)
    
    # If the process is still running after the timeout
    if process.is_alive():
        logger.warning(f"Retrosynthesis timed out after {timeout} seconds for SMILES: {smiles}")
        # Terminate the process
        process.terminate()
        process.join()
        return False
    
    # If the process finished within the timeout, get the result
    if not queue.empty():
        return queue.get()
    
    return False

def main(out_dir, checkpoint, pdbfile, resi_list, n_samples, sanitize,
         protein_file, receptor, program_choice="qvina", scoring_function="nnscore2",
         center=(114.817, 75.602, 82.416), box_size=(38, 70, 58),
         exhaustiveness=10, is_selfies=False, is_peptide=False, 
         top_n=5, max_variants=5):
    """
    Quick pipeline main function.
    
    Args:
        out_dir: Output directory
        checkpoint: Checkpoint file for ligand generation
        pdbfile: PDB file for ligand generation
        resi_list: Residue list for ligand generation
        n_samples: Number of samples for ligand generation
        sanitize: Whether to sanitize generated molecules
        protein_file: Protein file for docking
        receptor: Receptor file for docking
        program_choice: Docking program choice
        scoring_function: Scoring function for docking
        center: Center coordinates for docking box
        box_size: Box dimensions for docking
        exhaustiveness: Exhaustiveness for docking
        is_selfies: Whether to use SELFIES representation
        is_peptide: Whether the ligand is a peptide
        top_n: Number of top compounds to process (only used for final analysis)
        max_variants: Maximum number of variants per compound for retrosynthesis
    """
    # Set up output directories
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories
    ligand_gen_dir = out_dir / "ligand_generation"
    ligand_gen_dir.mkdir(exist_ok=True)
    
    retro_dir = out_dir / "retrosyn_results"
    retro_dir.mkdir(exist_ok=True)
    
    dock_dir = out_dir / "docking_results"
    dock_dir.mkdir(exist_ok=True)
    
    filter_dir = out_dir / "filter_results"
    filter_dir.mkdir(exist_ok=True)
    
    # Step 1: Ligand Generation
    logger.info("Running ligand generation...")
    base_name = "quick_pipeline"
    ligand_gen_out = ligand_gen_dir / f"{base_name}_mols_gen.sdf"
    
    lg_thread = run_ligand_generation(
        checkpoint=checkpoint,
        pdbfile=pdbfile,
        outfile=str(ligand_gen_out),
        resi_list=resi_list.split(),
        n_samples=n_samples,
        sanitize=sanitize,
        log_callback=logger.info
    )
    lg_thread.join()
    logger.info(f"Ligand generation complete. Output: {ligand_gen_out}")
    
    # Step 2: Extract SMILES from SDF
    compounds = extract_smiles_from_sdf(ligand_gen_out)
    if not compounds:
        logger.error("No valid compounds generated. Pipeline cannot continue.")
        return
    
    logger.info(f"Processing all {len(compounds)} valid generated compounds")
    
    # Add barcodes to each generated compound
    for idx, compound in enumerate(compounds):
        # Generate a unique barcode for each compound
        barcode = f"GEN1-{idx+1:04d}"
        compound["barcode"] = barcode
        compound["generation"] = "1"
    
    # Step 3: Run retrosynthesis on ALL compounds, not just top_n
    all_variants = []
    
    for idx, compound in enumerate(compounds):
        cid = compound["compound_id"]
        smiles = compound["smiles"]
        barcode = compound["barcode"]
        
        logger.info(f"Running retrosynthesis on compound {idx+1}/{len(compounds)}: {cid} (Barcode: {barcode})")
        
        # Create output file for retrosynthesis results
        retro_output = retro_dir / f"{cid}_retrosyn.csv"
        
        # Run retrosynthesis with a 5-minute timeout
        success = run_retrosynthesis_with_timeout(smiles, retro_output, timeout=300)
        
        if success:
            logger.info(f"Retrosynthesis completed for {cid}")
            
            # Extract variants from retrosynthesis results (top max_variants for each compound)
            variants = extract_variants_from_retrosynthesis(retro_output, max_variants=max_variants)
            
            # Add source compound information and barcodes to each variant
            for vidx, variant in enumerate(variants):
                # Barcode format: [Source Barcode]-V-[Variant Number]
                variant_barcode = f"{barcode}-V-{vidx+1:02d}"
                
                variant["source_compound"] = cid
                variant["source_smiles"] = smiles
                variant["barcode"] = variant_barcode
                variant["generation"] = "2"  # Generation 2 = retrosynthesis product
            
            all_variants.extend(variants)
        else:
            logger.warning(f"Retrosynthesis failed for {cid}")
    
    if not all_variants:
        logger.error("No variants generated from retrosynthesis. Pipeline cannot continue.")
        return
    
    logger.info(f"Total of {len(all_variants)} variants extracted from retrosynthesis")
    
    # Step 4: Apply MedChem filtering to variants
    filtered_variants = apply_medchem_filtering_to_variants(all_variants, filter_dir)
    
    if not filtered_variants:
        logger.warning("No variants passed MedChem filtering. Pipeline cannot continue.")
        return
    
    logger.info(f"After MedChem filtering, {len(filtered_variants)} variants remain")
    
    # Save filtered variants to SDF for reference
    filtered_sdf = filter_dir / "filtered_variants.sdf"
    smiles_to_sdf(filtered_variants, filtered_sdf)
    
    # Step 5: Redock filtered variants
    logger.info(f"Redocking {len(filtered_variants)} filtered variants...")
    
    # Prepare redocking parameters
    center_x, center_y, center_z = center
    size_x, size_y, size_z = box_size
    
    redock_params = (
        program_choice,
        scoring_function,
        center_x, center_y, center_z,
        size_x, size_y, size_z,
        exhaustiveness,
        is_selfies,
        is_peptide
    )
    
    # Create a directory for each variant's docking results
    redock_results = []
    for variant in filtered_variants:
        variant_id = variant["variant_id"]
        smiles = variant["smiles"]
        barcode = variant["barcode"]
        parent_id = variant.get("parent_id", "unknown")
        source_compound = variant.get("source_compound", "unknown")
        
        logger.info(f"Redocking variant: {variant_id} (Barcode: {barcode}, from {source_compound})")
        
        pose_out, rescored = redock_compound(
            variant_id,
            smiles,
            redock_params,
            receptor=receptor,
            log_callback=logger.info
        )
        
        redock_results.append({
            "variant_id": variant_id,
            "barcode": barcode,
            "parent_id": parent_id,
            "source_compound": source_compound,
            "smiles": smiles,
            "pose_pred_out": pose_out,
            "re_scored_values": rescored,
            "generation": "2"  # Generation 2 = retrosynthesis product
        })
        
        # Copy poses for this variant
        vfu_outputs_dir = Path(vfu_dir) / "outputs"
        variant_poses_dir = dock_dir / f"variant_{barcode}"
        variant_poses_dir.mkdir(exist_ok=True)
        
        # Copy all files from VFU outputs to the variant-specific directory
        for file_path in vfu_outputs_dir.glob("*"):
            if file_path.is_file():
                shutil.copy2(file_path, variant_poses_dir)
            elif file_path.is_dir():
                dest_dir = variant_poses_dir / file_path.name
                if dest_dir.exists():
                    shutil.rmtree(dest_dir)
                shutil.copytree(file_path, dest_dir)
        
        logger.info(f"Copied poses for variant {barcode} to {variant_poses_dir}")
    
    # Save results to CSV
    if redock_results:
        results_df = pd.DataFrame(redock_results)
        results_csv = dock_dir / "quick_pipeline_docking_results.csv"
        results_df.to_csv(results_csv, index=False)
        logger.info(f"Docking results saved to: {results_csv}")
        
        # Optionally select top_n compounds overall for final analysis
        if len(results_df) > top_n and 'pose_pred_out' in results_df.columns:
            try:
                # Extract best docking score and pose number for each compound
                def extract_best_pose_and_score(pose_pred_str):
                    try:
                        # Convert string representation to dictionary
                        import ast
                        
                        # Check if the pose_pred_str is already a dictionary
                        if isinstance(pose_pred_str, dict):
                            pose_dict = pose_pred_str
                        elif pose_pred_str is None:
                            return None, None
                        else:
                            # Try to safely evaluate the string representation
                            pose_dict = ast.literal_eval(pose_pred_str)
                        
                        # Find the pose with the best (lowest) score
                        best_score = float('inf')
                        best_pose = None
                        
                        for pose_name, pose_data in pose_dict.items():
                            # The pose_data is a list where the first element is a list of scores
                            # and the second element is the path to the pose file
                            if pose_data and isinstance(pose_data[0], list) and pose_data[0]:
                                score = pose_data[0][0]  # Get the first score from the first list
                                if score < best_score:
                                    best_score = score
                                    best_pose = pose_name
                        
                        if best_pose:
                            return best_score, best_pose.replace('.pdbqt', '')
                        else:
                            return None, None
                    
                    except Exception as e:
                        # Try to extract the best score directly from the string if parsing failed
                        try:
                            if isinstance(pose_pred_str, str):
                                # Look for patterns like [[-6.2, ...]] to extract the best score
                                import re
                                score_match = re.search(r'\[\[([-\d\.]+)', pose_pred_str)
                                if score_match:
                                    return float(score_match.group(1)), "unknown_pose"
                        except:
                            pass
                        
                        logger.warning(f"Error extracting pose data: {e}")
                        return None, None
                
                # Apply the extraction function to get best scores and poses
                results_df[['docking_score', 'best_pose']] = results_df['pose_pred_out'].apply(
                    lambda x: pd.Series(extract_best_pose_and_score(x))
                )
                
                # Filter out rows with None docking scores
                filtered_results = results_df.dropna(subset=['docking_score'])
                
                if len(filtered_results) > 0:
                    # Sort by docking score (lower is better) and select top_n unique compounds
                    top_results = filtered_results.sort_values('docking_score').head(top_n)
                    top_results_csv = dock_dir / f"quick_pipeline_top_{top_n}_results.csv"
                    top_results.to_csv(top_results_csv, index=False)
                    logger.info(f"Top {top_n} results saved to: {top_results_csv}")
                else:
                    logger.warning("No valid docking scores found after filtering")
            except Exception as e:
                logger.error(f"Error selecting top results: {e}")
    else:
        logger.warning("No docking results to save.")
    
    # Generate comprehensive tracking report
    logger.info("Generating compound tracking report...")
    generate_tracking_report(compounds, all_variants, redock_results, out_dir)
    
    logger.info("Quick pipeline completed successfully.")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Quick pipeline for drug discovery")
    
    # Required parameters
    parser.add_argument("--out_dir", type=str, help="Output directory", required=True)

    parser.add_argument("--checkpoint", type=str, default="DiffSBDD/checkpoints/crossdocked_fullatom_cond.ckpt",
                        help="Path to the checkpoint file (default: DiffSBDD/checkpoints/crossdocked_fullatom_cond.ckpt)")
    
    parser.add_argument("--pdbfile", type=str, default="input/NS5.pdb",
                        help="Path to target protein PDB file")
    
    parser.add_argument("--resi_list", type=str, default="A:719 A:770 A:841 A:856 A:887 A:888",
                        help="Residue identifiers (space-separated)")
    
    parser.add_argument("--protein_file", required=False, help="Protein file for docking", default="input/NS5_test.pdbqt")
    
    parser.add_argument("--receptor", required=False, help="Receptor file for docking", default="input/NS5_test.pdbqt")
    
    # Optional parameters with defaults
    parser.add_argument("--n_samples", type=int, default=10, help="Number of samples to generate")
    parser.add_argument("--sanitize", action="store_true", help="Sanitize generated molecules", default=True)
    parser.add_argument("--program_choice", default="qvina", help="Docking program choice")
    parser.add_argument("--scoring_function", default="nnscore2", help="Scoring function")
    parser.add_argument("--center", nargs=3, type=float, default=[114.817, 75.602, 82.416], 
                        help="Docking box center coordinates")
    parser.add_argument("--box_size", nargs=3, type=int, default=[38, 70, 58],
                        help="Docking box dimensions")
    parser.add_argument("--exhaustiveness", type=int, default=10, help="Docking exhaustiveness")
    parser.add_argument("--is_selfies", action="store_true", help="Use SELFIES representation")
    parser.add_argument("--is_peptide", action="store_true", help="Ligand is a peptide")
    parser.add_argument("--top_n", type=int, default=5, 
                        help="Number of top compounds to process")
    parser.add_argument("--max_variants", type=int, default=5,
                        help="Maximum number of variants per compound")
    
    args = parser.parse_args()
    
    main(
        args.out_dir,
        args.checkpoint,
        args.pdbfile,
        args.resi_list,
        args.n_samples,
        args.sanitize,
        args.protein_file,
        args.receptor,
        args.program_choice,
        args.scoring_function,
        args.center,
        args.box_size,
        args.exhaustiveness,
        args.is_selfies,
        args.is_peptide,
        args.top_n,
        args.max_variants
    )