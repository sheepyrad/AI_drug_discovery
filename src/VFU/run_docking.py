#!/usr/bin/env python3
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
