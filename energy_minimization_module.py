# energy_minimization_module.py
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolToMolFile
import os

logger = logging.getLogger("EnergyMinimization")
logger.propagate = True

def split_sdf_file(ligand_sdf: Path, output_dir: Path, base_name: str) -> list:
    output_dir.mkdir(parents=True, exist_ok=True)
    suppl = Chem.SDMolSupplier(str(ligand_sdf), sanitize=False)
    output_paths = []
    for idx, mol in enumerate(suppl):
        if mol is None:
            logger.warning(f"Skipping molecule {idx} as it could not be parsed.")
            continue
        try:
            ligand_file = output_dir / f"{base_name}_ligand_{idx}.sdf"
            MolToMolFile(mol, str(ligand_file))
            output_paths.append(ligand_file)
        except Exception as e:
            logger.warning(f"Error processing molecule {idx}: {e}. Skipping.")
    return output_paths

def optimize_ligand(protein_file, ligand_file, output_file, cache_dir, prep_only, name, platform_name="fastest", add_solvent=False):
    from energy_minimization import optimize_ligand_in_pocket  # assume this function exists
    return optimize_ligand_in_pocket(
        protein_file=protein_file,
        ligand_file=ligand_file,
        output_file=output_file,
        temp_dir=cache_dir,
        prep_only=prep_only,
        name=name,
        platform_name=platform_name,
        add_solvent=add_solvent
    )

def concatenate_sdf_files_sorted(input_folder, output_file):
    """
    Concatenates all SDF files in the input_folder into a single SDF file.
    The files are sorted in ascending order based on the numeric value extracted
    from their filename using the pattern "_ligand_(\\d+)_minimized.sdf".
    """
    import glob, re
    input_folder = str(input_folder)
    sdf_files = glob.glob(os.path.join(input_folder, "*.sdf"))
    if not sdf_files:
        logger.info("No SDF files found in the specified folder.")
        return None

    def extract_number(fname):
        m = re.search(r"_ligand_(\d+)_minimized\.sdf$", fname)
        return int(m.group(1)) if m else float('inf')

    sdf_files_sorted = sorted(sdf_files, key=lambda f: extract_number(os.path.basename(f)))
    logger.info("Concatenation order:")
    for f in sdf_files_sorted:
        logger.info(f"{os.path.basename(f)} (extracted number: {extract_number(os.path.basename(f))})")

    writer = Chem.SDWriter(str(output_file))
    for sdf_file in sdf_files_sorted:
        logger.info(f"Processing file: {sdf_file}")
        supplier = Chem.SDMolSupplier(sdf_file)
        for mol in supplier:
            if mol is not None:
                writer.write(mol)
    writer.close()
    logger.info(f"Concatenated SDF file saved to: {output_file}")
    return output_file
