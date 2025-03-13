# medchem_filter.py
import os, re
import pandas as pd
import numpy as np
from rdkit import Chem
import medchem as mc
import logging
from pathlib import Path
import gc
import time

logger = logging.getLogger("MedChemFilter")
logger.propagate = True

def extract_number_from_string(s):
    numbers = re.findall(r'\d+', s)
    return int(numbers[-1]) if numbers else 0

def load_compounds(sdf_file):
    supplier = Chem.SDMolSupplier(str(sdf_file))
    if supplier is None or len(supplier) == 0:
        raise ValueError("No molecules were loaded from the SDF file. Check its content.")
    compounds = []
    for i, mol in enumerate(supplier):
        if mol is None:
            continue
        compound_id = mol.GetProp("_Name") if mol.HasProp("_Name") else f"compound_{i}"
        file_num = extract_number_from_string(compound_id)
        smiles = Chem.MolToSmiles(mol)
        compounds.append({
            "mol": mol,
            "smiles": smiles,
            "compound_id": compound_id,
            "mol_index": file_num
        })
    if not compounds:
        raise ValueError("No compounds were loaded from the SDF file. Check its contents.")
    df = pd.DataFrame(compounds)
    df["file_index"] = df["mol_index"].astype(str)
    return df

# Add a function to process molecules in batches to avoid file descriptor exhaustion
def process_in_batches(molecules, func, batch_size=10, **kwargs):
    """
    Process molecules in small batches to avoid file descriptor exhaustion.
    
    Args:
        molecules: List of RDKit molecules
        func: Function to apply to each batch
        batch_size: Size of each batch
        **kwargs: Additional arguments to pass to func
        
    Returns:
        List of results
    """
    results = []
    total = len(molecules)
    
    for i in range(0, total, batch_size):
        batch = molecules[i:min(i+batch_size, total)]
        logger.info(f"Processing batch {i//batch_size + 1}/{(total+batch_size-1)//batch_size} ({len(batch)} molecules)")
        
        # Apply the function to the batch
        batch_results = func(mols=batch, **kwargs)
        results.extend(batch_results)
        
        # Force garbage collection and sleep briefly to allow file descriptors to be released
        gc.collect()
        time.sleep(1)
    
    return results

def apply_filters(df, threshold=12):
    # Apply MedChem rule-based and functional filters.
    logger.info("Applying basic rule-based filters...")
    df["rule_of_five"]                = df["smiles"].apply(mc.rules.basic_rules.rule_of_five)
    df["rule_of_ghose"]               = df["smiles"].apply(mc.rules.basic_rules.rule_of_ghose)
    df["rule_of_veber"]               = df["smiles"].apply(mc.rules.basic_rules.rule_of_veber)
    df["rule_of_reos"]                = df["smiles"].apply(mc.rules.basic_rules.rule_of_reos)
    df["rule_of_chemaxon_druglikeness"] = df["smiles"].apply(mc.rules.basic_rules.rule_of_chemaxon_druglikeness)
    df["rule_of_egan"]                = df["smiles"].apply(mc.rules.basic_rules.rule_of_egan)
    df["rule_of_pfizer_3_75"]         = df["smiles"].apply(mc.rules.basic_rules.rule_of_pfizer_3_75)
    df["rule_of_gsk_4_400"]           = df["smiles"].apply(mc.rules.basic_rules.rule_of_gsk_4_400)
    df["rule_of_oprea"]               = df["smiles"].apply(mc.rules.basic_rules.rule_of_oprea)
    df["rule_of_xu"]                  = df["smiles"].apply(mc.rules.basic_rules.rule_of_xu)
    df["rule_of_zinc"]                = df["smiles"].apply(mc.rules.basic_rules.rule_of_zinc)
    df["rule_of_leadlike_soft"]       = df["smiles"].apply(mc.rules.basic_rules.rule_of_leadlike_soft)
    df["rule_of_druglike_soft"]       = df["smiles"].apply(mc.rules.basic_rules.rule_of_druglike_soft)
    df["rule_of_generative_design"]   = df["smiles"].apply(mc.rules.basic_rules.rule_of_generative_design)
    df["rule_of_generative_design_strict"] = df["smiles"].apply(mc.rules.basic_rules.rule_of_generative_design_strict)
    
    # Limit the number of parallel jobs to avoid exhausting file descriptors
    n_jobs = 1  # Use only 1 job to avoid "Too many open files" error
    
    # Get the list of molecules
    molecules = df["mol"].tolist()
    
    # Process each filter in small batches with garbage collection in between
    logger.info("Applying BMS alerts filter in batches...")
    try:
        df["alerts_BMS"] = process_in_batches(
            molecules, 
            mc.functional.alert_filter, 
            batch_size=5,  # Very small batch size to avoid file descriptor issues
            alerts=["BMS"], 
            n_jobs=n_jobs, 
            progress=False, 
            return_idx=False
        )
        gc.collect()
    except Exception as e:
        logger.error(f"Error in BMS alerts filter: {e}")
        df["alerts_BMS"] = [False] * len(df)
    
    logger.info("Applying PAINS alerts filter in batches...")
    try:
        df["alerts_PAINS"] = process_in_batches(
            molecules, 
            mc.functional.alert_filter, 
            batch_size=5,
            alerts=["PAINS"], 
            n_jobs=n_jobs, 
            progress=False, 
            return_idx=False
        )
        gc.collect()
    except Exception as e:
        logger.error(f"Error in PAINS alerts filter: {e}")
        df["alerts_PAINS"] = [False] * len(df)
    
    logger.info("Applying SureChEMBL alerts filter in batches...")
    try:
        df["alerts_SureChEMBL"] = process_in_batches(
            molecules, 
            mc.functional.alert_filter, 
            batch_size=5,
            alerts=["SureChEMBL"], 
            n_jobs=n_jobs, 
            progress=False, 
            return_idx=False
        )
        gc.collect()
    except Exception as e:
        logger.error(f"Error in SureChEMBL alerts filter: {e}")
        df["alerts_SureChEMBL"] = [False] * len(df)
    
    logger.info("Applying NIBR filter in batches...")
    try:
        df["filters_NIBR"] = process_in_batches(
            molecules, 
            mc.functional.nibr_filter, 
            batch_size=5,
            n_jobs=n_jobs, 
            progress=False, 
            return_idx=False
        )
        gc.collect()
    except Exception as e:
        logger.error(f"Error in NIBR filter: {e}")
        df["filters_NIBR"] = [False] * len(df)
    
    logger.info("Applying complexity filter in batches...")
    try:
        df["filter_complexity"] = process_in_batches(
            molecules, 
            mc.functional.complexity_filter, 
            batch_size=5,
            complexity_metric="bertz", 
            threshold_stats_file="zinc_15_available", 
            n_jobs=n_jobs, 
            progress=False, 
            return_idx=False
        )
        gc.collect()
    except Exception as e:
        logger.error(f"Error in complexity filter: {e}")
        df["filter_complexity"] = [False] * len(df)
    
    logger.info("Applying bredt filter in batches...")
    try:
        df["filter_bredt"] = process_in_batches(
            molecules, 
            mc.functional.bredt_filter, 
            batch_size=5,
            n_jobs=n_jobs, 
            progress=False, 
            return_idx=False
        )
        gc.collect()
    except Exception as e:
        logger.error(f"Error in bredt filter: {e}")
        df["filter_bredt"] = [False] * len(df)
    
    logger.info("Applying molecular graph filter in batches...")
    try:
        df["filter_molecular_graph"] = process_in_batches(
            molecules, 
            mc.functional.molecular_graph_filter, 
            batch_size=5,
            max_severity=5, 
            n_jobs=n_jobs, 
            progress=False, 
            return_idx=False
        )
        gc.collect()
    except Exception as e:
        logger.error(f"Error in molecular graph filter: {e}")
        df["filter_molecular_graph"] = [False] * len(df)
    
    logger.info("Applying lilly demerit filter in batches...")
    try:
        df["filter_lilly_demerit"] = process_in_batches(
            molecules, 
            mc.functional.lilly_demerit_filter, 
            batch_size=5,
            n_jobs=n_jobs, 
            progress=False, 
            return_idx=False
        )
        gc.collect()
    except Exception as e:
        logger.error(f"Error in lilly demerit filter: {e}")
        df["filter_lilly_demerit"] = [False] * len(df)
    
    filter_columns = [
        "rule_of_five",
        "rule_of_ghose",
        "rule_of_veber",
        "rule_of_reos",
        "rule_of_chemaxon_druglikeness",
        "rule_of_egan",
        "rule_of_pfizer_3_75",
        "rule_of_gsk_4_400",
        "rule_of_oprea",
        "rule_of_xu",
        "rule_of_zinc",
        "rule_of_leadlike_soft",
        "rule_of_druglike_soft",
        "rule_of_generative_design",
        "rule_of_generative_design_strict",
        "alerts_BMS",
        "alerts_PAINS",
        "alerts_SureChEMBL",
        "filters_NIBR",
        "filter_complexity",
        "filter_bredt",
        "filter_molecular_graph",
        "filter_lilly_demerit",
    ]
    for col in filter_columns:
        df[col] = df[col].apply(lambda x: 1 if x else 0)
    df["n_filters_pass"] = df[filter_columns].sum(axis=1)
    df["gen_design_score"] = df["rule_of_generative_design"] + df["rule_of_generative_design_strict"]
    
    return df

def generative_filter(sdf_file, output_folder=None):
    """
    Loads compounds from an SDF file, applies all MedChem filters, and then
    returns only those compounds that pass at least one of the generative design tests.
    Also, outputs two CSV files:
      1. All compounds with computed filter values.
      2. The filtered compounds (passing generative criteria).
      
    Args:
        sdf_file: Path to the SDF file containing compounds
        output_folder: Directory to save results (if None, will use parent directory of sdf_file)
        
    Returns:
        DataFrame with filtered compounds
    """
    try:
        logger.info(f"Loading compounds from {sdf_file}...")
        df = load_compounds(sdf_file)
        logger.info(f"Loaded {len(df)} compounds. Applying filters...")
        df = apply_filters(df)
        gc.collect()
        
        # Use the parent directory of the sdf_file if output_folder is not specified
        if output_folder is None:
            output_folder = Path(sdf_file).parent / "medchem_results"
        
        # Ensure the output folder exists.
        output_folder = Path(output_folder)
        output_folder.mkdir(exist_ok=True, parents=True)
        
        # Save unfiltered compounds.
        pre_filter_file = output_folder / "all_compounds.csv"
        df.to_csv(pre_filter_file, index=False)
        logger.info(f"Saved all compounds (with filter metrics computed) to: {pre_filter_file}")
        
        # Apply the generative filter.
        df_filtered = df[
            (df["rule_of_generative_design"] == 1) | 
            (df["rule_of_generative_design_strict"] == 1)
        ]
        filtered_file = output_folder / "generative_filter_results.csv"
        df_filtered.to_csv(filtered_file, index=False)
        logger.info(f"Generative filtering: {len(df_filtered)} compounds remain. Results saved to: {filtered_file}")
        
        return df_filtered
    except Exception as e:
        logger.error(f"Error in generative_filter: {e}")
        # Return an empty DataFrame if there's an error
        return pd.DataFrame()

def filter_compounds(sdf_file):
    # Simply call generative_filter with the SDF file path.
    return generative_filter(sdf_file)

