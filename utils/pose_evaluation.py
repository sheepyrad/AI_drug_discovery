# pose_evaluation.py
import pandas as pd
from posebusters import PoseBusters
import os
from pathlib import Path

def run_posebuster(sdf_file, config="mol", output_folder=None):
    """
    Run PoseBusters on the input SDF file.
    
    Args:
        sdf_file: Path to the SDF file to analyze
        config: Configuration to use for PoseBusters
        output_folder: Directory to save results (if None, will use parent directory of sdf_file)
        
    Returns:
        DataFrame with PoseBusters results
    """
    buster = PoseBusters(config=config)
    df = buster.bust(sdf_file, None, None)
    df['all_true'] = df.all(axis=1)
    
    # Use the parent directory of the sdf_file if output_folder is not specified
    if output_folder is None:
        output_folder = Path(sdf_file).parent / "posebuster_results"
    
    # Save the DataFrame to a CSV file
    output_folder = Path(output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)
    output_file = output_folder / "posebuster_results.csv"
    df.to_csv(output_file, index=False)
    
    return df

def extract_valid_ligands(df_ligands, minimized_ligands_folder, output_sdf_file):
    from rdkit import Chem
    import os
    writer = Chem.SDWriter(str(output_sdf_file))
    valid_ligand_names = sorted(df_ligands[df_ligands["all_true"]].index.get_level_values("molecule"),
                                  key=lambda name: int(name.split("_")[-1]))
    count = 0
    for ligand_name in valid_ligand_names:
        sdf_file = minimized_ligands_folder / f"{ligand_name}_minimized.sdf"
        if sdf_file.exists():
            supplier = Chem.SDMolSupplier(str(sdf_file))
            if supplier and supplier[0] is not None:
                mol = supplier[0]
                mol_name = mol.GetProp("_Name").strip() if mol.HasProp("_Name") else ""
                if mol_name == ligand_name:
                    writer.write(mol)
                    count += 1
    writer.close()
    return output_sdf_file, count
