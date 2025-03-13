# redocking.py
import os
import sys
from pathlib import Path
import shutil
import logging

# Determine the absolute path to the VFU folder relative to this file
vfu_dir = Path(__file__).parent / "VFU"
vfu_config_dir = vfu_dir / "config"
input_dir = Path(__file__).parent / "input"

def redock_compound(compound_id, smiles, redock_params, receptor, log_callback=None):
    """
    Redock a compound using VF Unity.
    
    Args:
        compound_id: ID of the compound
        smiles: SMILES string of the compound
        redock_params: Tuple of redocking parameters
        receptor: Filename of the receptor file (should be in ./input directory)
        log_callback: Function to call for logging
    
    Returns:
        Tuple of (pose_pred_out, re_scored_values)
    """
    (program_choice, scoring_function, center_x, center_y, center_z,
     size_x, size_y, size_z, exhaustiveness, is_selfies, is_peptide) = redock_params

    # Get just the filename from the receptor path
    receptor_filename = os.path.basename(receptor)
    
    # Expected location in input directory
    receptor_in_input = input_dir / receptor_filename
    
    # Final path where the receptor should be in VFU/config
    receptor_in_vfu_config = vfu_config_dir / receptor_filename
    
    if log_callback:
        log_callback(f"Redocking {compound_id} with SMILES: {smiles}")
        log_callback(f"Looking for receptor file: {receptor_filename} in input directory")
    
    # Check if receptor exists in input directory
    if not receptor_in_input.exists():
        error_msg = f"Receptor file '{receptor_filename}' not found in input directory. Please place it in the ./input directory."
        if log_callback:
            log_callback(error_msg)
        raise FileNotFoundError(error_msg)
    
    # Store current directory
    current_dir = os.getcwd()
    
    # Ensure VFU config directory exists
    if not vfu_config_dir.exists():
        if log_callback:
            log_callback(f"Creating VFU config directory: {vfu_config_dir}")
        vfu_config_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy receptor from input to VFU/config if needed
    if not receptor_in_vfu_config.exists():
        if log_callback:
            log_callback(f"Copying receptor from {receptor_in_input} to {receptor_in_vfu_config}")
        shutil.copy2(receptor_in_input, receptor_in_vfu_config)
    
    try:
        # Change to VFU directory
        os.chdir(vfu_dir)
        
        # Add VFU directory to path temporarily
        if str(vfu_dir) not in sys.path:
            sys.path.insert(0, str(vfu_dir))
        
        # Import here after changing directory and path
        from run_vf_unity import main
        
        # Use relative path for the receptor after changing directories
        receptor_rel_path = f"./config/{receptor_filename}"
        
        if log_callback:
            log_callback(f"Using receptor path for VFU: {receptor_rel_path}")
        
        pose_pred_out, re_scored_values = main(
            program_choice,
            scoring_function,
            center_x,
            center_y,
            center_z,
            size_x,
            size_y,
            size_z,
            exhaustiveness,
            smiles,
            is_selfies,
            is_peptide,
            receptor_rel_path  # Use relative path here
        )
        
        return pose_pred_out, re_scored_values
        
    except Exception as e:
        if log_callback:
            log_callback(f"Error redocking {compound_id}: {e}")
        return None, None
    
    finally:
        # Change back to original directory
        os.chdir(current_dir)
        
        # Remove VFU directory from path if we added it
        if str(vfu_dir) in sys.path:
            sys.path.remove(str(vfu_dir))