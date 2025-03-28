# redocking.py
import os
import sys
from pathlib import Path
import shutil
import logging
import traceback

# Determine the absolute path to the VFU folder relative to this file
vfu_dir = Path(__file__).parent.parent / "src" / "VFU"
vfu_config_dir = vfu_dir / "config"
vfu_inputs_dir = vfu_dir / "inputs"
vfu_outputs_dir = vfu_dir / "outputs"
input_dir = Path(__file__).parent.parent / "input"

# Get logger for this module
logger = logging.getLogger(__name__)

def redock_compound(compound_id, smiles, redock_params, receptor=None, log_callback=print):
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
    if log_callback is None:
        log_callback = print
    
    log_callback(f"Starting redocking process for {compound_id}")
    
    # Check if VFU directory exists
    if not vfu_dir.exists():
        error_msg = f"VFU directory not found at {vfu_dir}. Please ensure it exists."
        log_callback(error_msg)
        return None, None
        
    # Ensure VFU directories exist
    vfu_config_dir.mkdir(parents=True, exist_ok=True)
    vfu_inputs_dir.mkdir(parents=True, exist_ok=True)
    vfu_outputs_dir.mkdir(parents=True, exist_ok=True)
    
    log_callback(f"VFU directories initialized: config={vfu_config_dir.exists()}, inputs={vfu_inputs_dir.exists()}, outputs={vfu_outputs_dir.exists()}")
        
    (program_choice, scoring_function, center_x, center_y, center_z,
     size_x, size_y, size_z, exhaustiveness, is_selfies, is_peptide) = redock_params

    # Get just the filename from the receptor path if provided
    receptor_filename = None
    if receptor:
        receptor_filename = os.path.basename(receptor)
        
        # Expected location in input directory
        receptor_in_input = input_dir / receptor_filename
        
        # Final path where the receptor should be in VFU/config
        receptor_in_vfu_config = vfu_config_dir / receptor_filename
        
        log_callback(f"Redocking {compound_id} with SMILES: {smiles}")
        log_callback(f"Looking for receptor file: {receptor_filename} in input directory")
    
        # Check if receptor exists in input directory
        if not receptor_in_input.exists():
            error_msg = f"Receptor file '{receptor_filename}' not found in input directory at {receptor_in_input}. Please place it there."
            log_callback(error_msg)
            return None, None
    
        # Copy receptor from input to VFU/config if needed
        if not receptor_in_vfu_config.exists():
            log_callback(f"Copying receptor from {receptor_in_input} to {receptor_in_vfu_config}")
            try:
                shutil.copy2(receptor_in_input, receptor_in_vfu_config)
                log_callback(f"Receptor copied successfully. File exists: {receptor_in_vfu_config.exists()}")
            except Exception as e:
                log_callback(f"Error copying receptor file: {e}")
                return None, None

    # Instead of changing directory, we'll modify the config.txt file in the VFU directory
    # and run a custom function that mimics the behavior of run_vf_unity.main
    
    # Create a custom config file
    config_file = vfu_dir / "config.txt"
    log_callback(f"Creating custom config file at {config_file}")
    
    try:
        with open(config_file, 'w') as f:
            f.write(f"program_choice={program_choice}\n")
            f.write(f"center_x={center_x}\n")
            f.write(f"center_y={center_y}\n")
            f.write(f"center_z={center_z}\n")
            f.write(f"size_x={size_x}\n")
            f.write(f"size_y={size_y}\n")
            f.write(f"size_z={size_z}\n")
            f.write(f"exhaustiveness={exhaustiveness}\n")
            f.write(f"smi={smiles}\n")
            f.write(f"is_selfies={str(is_selfies)}\n")
            f.write(f"is_peptide={str(is_peptide)}\n")
            
            # Set receptor path relative to VFU directory
            if receptor_filename:
                f.write(f"receptor=./config/{receptor_filename}\n")
            else:
                f.write("receptor=\n")
                
        log_callback(f"Config file created successfully")
    except Exception as e:
        log_callback(f"Error creating config file: {e}")
        return None, None
    
    # Create a wrapper script in the VFU directory
    wrapper_script = vfu_dir / "run_docking.py"
    log_callback(f"Creating wrapper script at {wrapper_script}")
    
    try:
        with open(wrapper_script, 'w') as f:
            f.write("""#!/usr/bin/env python3
import os
import sys
import json
from run_vf_unity import main, read_config_file

if __name__ == "__main__":
    # Read configuration parameters
    program_choice, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, is_selfies, is_peptide, receptor = read_config_file()
    
    # Handle scoring function
    scoring_function = ""
    if '+' in program_choice:
        program_choice, scoring_function = program_choice.split('+')[0], program_choice.split('+')[1]
    
    # Run the main function
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
        smi,
        is_selfies,
        is_peptide,
        receptor
    )
    
    # Output results in a format that can be captured
    result = {
        "pose_pred_out": pose_pred_out,
        "re_scored_values": re_scored_values
    }
    print("RESULT_JSON_START")
    print(json.dumps(result))
    print("RESULT_JSON_END")
""")
        log_callback(f"Wrapper script created successfully")
    except Exception as e:
        log_callback(f"Error creating wrapper script: {e}")
        return None, None
    
    # Run the wrapper script from the VFU directory without changing our working directory
    import subprocess
    import json
    
    try:
        log_callback(f"Running docking script...")
        # Run the script with Python - use the same Python executable that's running this code
        python_exec = sys.executable
        cmd = [python_exec, str(wrapper_script)]
        
        # Set the working directory in the subprocess call
        result = subprocess.run(
            cmd, 
            cwd=str(vfu_dir),
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            log_callback(f"Docking failed with error code {result.returncode}")
            log_callback(f"STDOUT: {result.stdout}")
            log_callback(f"STDERR: {result.stderr}")
            return None, None
        
        # Parse the output to extract the JSON results
        stdout = result.stdout
        log_callback(f"Docking completed, parsing results...")
        
        # Extract JSON results
        start_marker = "RESULT_JSON_START"
        end_marker = "RESULT_JSON_END"
        
        if start_marker in stdout and end_marker in stdout:
            json_start = stdout.find(start_marker) + len(start_marker)
            json_end = stdout.find(end_marker)
            json_str = stdout[json_start:json_end].strip()
            
            try:
                docking_results = json.loads(json_str)
                pose_pred_out = docking_results.get("pose_pred_out")
                re_scored_values = docking_results.get("re_scored_values")
                
                log_callback(f"Docking results parsed successfully")
                return pose_pred_out, re_scored_values
            except json.JSONDecodeError as e:
                log_callback(f"Error parsing JSON results: {e}")
                log_callback(f"Raw JSON string: {json_str}")
                return None, None
        else:
            log_callback(f"Couldn't find result markers in output")
            log_callback(f"STDOUT: {stdout}")
            return None, None
            
    except Exception as e:
        log_callback(f"Error executing docking script: {str(e)}")
        log_callback(f"Traceback: {traceback.format_exc()}")
        return None, None