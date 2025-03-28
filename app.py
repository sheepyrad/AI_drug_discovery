import streamlit as st
import os
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol

# Set page config
st.set_page_config(
    page_title="Drug Discovery Pipeline",
    page_icon="💊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Add custom CSS
st.markdown("""
    <style>
    .main {
        padding: 0rem 1rem;
    }
    .stApp {
        max-width: 100%;
        margin: 0 auto;
    }
    .block-container {
        max-width: 100%;
        padding-left: 1rem;
        padding-right: 1rem;
    }
    .st-emotion-cache-1v0mbdj {
        width: 100%;
    }
    </style>
""", unsafe_allow_html=True)

# Title and description
st.title("🧬 Drug Discovery Pipeline")
st.markdown("""
    ### Welcome to the Drug Discovery Pipeline Application
    
    This application provides an interactive interface for running and monitoring our 
    drug discovery pipeline. The pipeline includes multiple stages:
    
    1. **Ligand Generation** 🔬
    2. **Retrosynthesis Analysis** 🧪
    3. **MedChem Filtering** ⚗️
    4. **Molecular Docking** 🎯
    
    Navigate through the pages using the sidebar to:
    - Configure pipeline parameters
    - Execute the pipeline
    - View and analyze results
    - Visualize existing results from previous runs
""")


# We don't need custom navigation as Streamlit multipage apps have built-in navigation

# Quick start guide
st.header("Quick Start Guide")
with st.expander("How to use this application"):
    st.markdown("""
        1. **Configuration Page**
           - Upload required files (PDB, checkpoint)
           - Set pipeline parameters
           - Validate inputs
        
        2. **Execution Page**
           - Start pipeline execution
           - Monitor progress
           - View real-time logs
        
        3. **Results Page**
           - Visualize compounds
           - Analyze docking results
           - Export data
           
        4. **Visualize Results Page**
           - Load results from existing output directories
           - View and analyze previous run results
           - Download pose files and structures
    """)

# System status
st.header("System Status")
col1, col2, col3 = st.columns(3)

with col1:
    st.metric(
        label="GPU Status",
        value="Available" if os.environ.get("CUDA_VISIBLE_DEVICES") else "Not Available"
    )

with col2:
    st.metric(
        label="Memory Usage",
        value=f"{os.popen('free -h').readlines()[1].split()[2]}"
    )

with col3:
    st.metric(
        label="Storage",
        value=f"{os.popen('df -h /').readlines()[1].split()[3]}"
    )

# Recent runs
st.header("Recent Pipeline Runs")
if "recent_runs" not in st.session_state:
    st.session_state.recent_runs = []

if st.session_state.recent_runs:
    runs_df = pd.DataFrame(st.session_state.recent_runs)
    st.dataframe(runs_df)
else:
    st.info("No recent pipeline runs found. Start a new run from the Configuration page or visualize existing results from the Visualize Results page.")

# Footer
st.markdown("---")
st.markdown("""
    <div style='text-align: center'>
        <p>Drug Discovery Pipeline Application | Version 1.0.0</p>
        <p>For support or issues, please visit our <a href='https://github.com/yourusername/drug_pipeline_app'>GitHub repository</a></p>
    </div>
""", unsafe_allow_html=True)

# Initialize session state variables if they don't exist
if "pipeline_config" not in st.session_state:
    st.session_state.pipeline_config = {
        "out_dir": None,
        "checkpoint": None,
        "pdbfile": None,
        "resi_list": None,
        "n_samples": 200,
        "sanitize": True,
        "protein_file": None,
        "receptor": None,
        "program_choice": "qvina",
        "scoring_function": "nnscore2",
        "center": [114.817, 75.602, 82.416],
        "box_size": [38, 70, 58],
        "exhaustiveness": 32,
        "is_selfies": False,
        "is_peptide": False,
        "top_n": 5,
        "max_variants": 5,
        "num_rounds": 1
    }

if "pipeline_status" not in st.session_state:
    st.session_state.pipeline_status = {
        "running": False,
        "current_round": 0,
        "total_rounds": 0,
        "current_stage": None,
        "progress": 0,
        "log_messages": []
    }

if "results" not in st.session_state:
    st.session_state.results = {
        "compounds": [],
        "variants": [],
        "docking_results": [],
        "tracking_report": None
    } 