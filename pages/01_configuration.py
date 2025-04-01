import streamlit as st
import os
from pathlib import Path
import tempfile
import py3Dmol
import re
import streamlit.components.v1 as components
import logging

# Set up logging
logger = logging.getLogger(__name__)

st.set_page_config(
    page_title="Pipeline Configuration",
    page_icon="⚙️",
    layout="wide"
)

# Add custom CSS for full-width layout
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

# Initialize session state variables
if "pipeline_config" not in st.session_state:
    st.session_state.pipeline_config = None

if "stop_pipeline" not in st.session_state:
    st.session_state.stop_pipeline = False

def parse_residue_list(resi_list):
    """Parse residue list string into a list of chain and residue numbers"""
    residues = []
    for res in resi_list.split():
        if ':' in res:
            chain, num = res.split(':')
            residues.append((chain, int(num)))
    return residues

def add_box(view, center=[0,0,0], dimensions=[10,10,10], color='blue', opacity=0.4, add_wireframe=True):
    """Adds a box to a py3Dmol view object.
    
    Parameters
    ----------
    view : py3Dmol.view
        The py3Dmol view object to add the box to
    center : list of float, default [0,0,0]
        The x,y,z coordinates of the box center in Angstroms
    dimensions : list of float, default [10,10,10]
        The width, height, depth of the box in Angstroms
    color : str, default 'blue'
        The color of the box
    opacity : float, default 0.4
        The opacity of the box (0.0 to 1.0)
    add_wireframe : bool, default True
        Whether to add a wireframe outline of the box
    """
    # Add transparent box
    view.addBox({
        'center': {'x': center[0], 'y': center[1], 'z': center[2]},
        'dimensions': {'w': dimensions[0], 'h': dimensions[1], 'd': dimensions[2]},
        'color': color,
        'opacity': opacity
    })
    
    # Add wireframe box if requested
    if add_wireframe:
        view.addBox({
            'center': {'x': center[0], 'y': center[1], 'z': center[2]},
            'dimensions': {'w': dimensions[0], 'h': dimensions[1], 'd': dimensions[2]},
            'color': color,
            'wireframe': True
        })

def visualize_protein_residues(pdb_content, selected_residues, center=None, box_size=None):
    """Create a py3Dmol visualization of the protein with highlighted residues and docking box"""
    view = py3Dmol.view(width=800, height=600)
    view.addModel(pdb_content, "pdb")
    
    # Set up the basic protein visualization
    view.setStyle({'cartoon': {'color': 'lightgray'}})
    
    # Highlight selected residues
    for chain, resnum in selected_residues:
        view.addStyle({
            'chain': chain,
            'resi': resnum
        }, {
            'cartoon': {'color': 'red'},
            'stick': {'color': 'red'},
            'labels': {'fontColor': 'black', 'showResname': True}
        })
    
    # Add docking box if parameters are provided
    if center and box_size:
        add_box(view, center=center, dimensions=box_size)
    
    view.zoomTo()
    
    # Get the HTML representation
    html = view._make_html()
    
    # Return both the view object and HTML
    return view, html

st.title("⚙️ Pipeline Configuration")
st.markdown("""
    Configure the parameters for your drug discovery pipeline run.
    Required fields are marked with an asterisk (*).
""")

# Create form for configuration
with st.form("pipeline_config_form"):
    st.header("File Uploads")
    
    # File uploads
    pdb_file = st.file_uploader(
        "PDB File *",
        type=["pdb"],
        help="Upload the target protein PDB file"
    )
    
    receptor_file = st.file_uploader(
        "Receptor File",
        type=["pdbqt"],
        help="Upload the receptor file for docking (optional)"
    )
    
    if not receptor_file:
        st.warning("⚠️ No receptor file uploaded. Docking steps may be skipped or fail. Upload a PDBQT receptor file if you want to perform docking.")
    
    st.header("Basic Parameters")
    
    col1, col2 = st.columns(2)
    
    with col1:
        resi_list = st.text_input(
            "Residue List *",
            value="A:719 A:770 A:841 A:856 A:887 A:888",
            help="Space-separated residue identifiers"
        )
        
        n_samples = st.number_input(
            "Number of Samples",
            min_value=1,
            max_value=500,
            value=5,
            help="Number of compounds to generate"
        )
        
        sanitize = st.checkbox(
            "Sanitize Generated Molecules",
            value=True,
            help="Apply sanitization to generated molecules"
        )
    
    with col2:
        program_choice = st.selectbox(
            "Docking Program",
            options=["qvina"],
            index=0,
            help="Select the docking program to use"
        )
        
        scoring_function = st.selectbox(
            "Scoring Function",
            options=["nnscore2"],
            index=0,
            help="Select the scoring function for docking"
        )
    
    st.header("Advanced Parameters")
    
    # Docking box parameters
    st.subheader("Docking Box Configuration")
    col3, col4 = st.columns(2)
    
    with col3:
        st.write("Center Coordinates")
        center_x = st.number_input("Center X", value=114.817, format="%.3f")
        center_y = st.number_input("Center Y", value=75.602, format="%.3f")
        center_z = st.number_input("Center Z", value=82.416, format="%.3f")
    
    with col4:
        st.write("Box Dimensions")
        size_x = st.number_input("Size X", value=38, min_value=1)
        size_y = st.number_input("Size Y", value=70, min_value=1)
        size_z = st.number_input("Size Z", value=58, min_value=1)
    
    # Add visualization section if PDB is uploaded
    if pdb_file is not None:
        st.write("### Protein Structure Visualization")
        st.write("Selected residues are highlighted in red. Docking box shown in blue.")
        
        # Read PDB content
        pdb_content = pdb_file.getvalue().decode('utf-8')
        
        # Parse residue list and create visualization
        selected_residues = parse_residue_list(resi_list)
        
        # Get current box parameters
        center = [center_x, center_y, center_z]
        box_size = [size_x, size_y, size_z]
        
        # Create visualization with box
        view, html = visualize_protein_residues(pdb_content, selected_residues, center, box_size)
        
        # Display the visualization using HTML component
        components.html(html, height=600, width=800)
        
        # Add some helpful instructions
        st.caption("""
        **Controls:**
        - Rotate: Click and drag
        - Zoom: Scroll wheel
        - Pan: Right click and drag
        - Reset view: Double click
        
        **Visualization Guide:**
        - Protein backbone: Light gray cartoon
        - Selected residues: Red cartoon and sticks with labels
        - Docking box: Blue transparent box with wireframe
        """)
    
    # Other parameters
    col5, col6 = st.columns(2)
    
    with col5:
        exhaustiveness = st.number_input(
            "Exhaustiveness",
            min_value=1,
            max_value=100,
            value=32,
            help="Docking exhaustiveness parameter"
        )
        
        is_selfies = st.checkbox(
            "Use SELFIES Representation",
            value=False,
            help="Use SELFIES molecular representation"
        )
    
    with col6:
        is_peptide = st.checkbox(
            "Ligand is Peptide",
            value=False,
            help="Check if the ligand is a peptide"
        )
        
        num_rounds = st.number_input(
            "Number of Rounds",
            min_value=1,
            max_value=100000,
            value=1,
            help="Number of pipeline rounds to run"
        )
    
    # Output configuration
    st.header("Output Configuration")
    
    col7, col8 = st.columns(2)
    
    with col7:
        output_dir_name = st.text_input(
            "Output Directory Name",
            value="pipeline_output",
            help="Name of the output directory where results will be stored"
        )
        
        top_n = st.number_input(
            "Top N Compounds",
            min_value=1,
            max_value=100,
            value=5,
            help="Number of top compounds to process"
        )
    
    with col8:
        max_variants = st.number_input(
            "Maximum Variants per Compound",
            min_value=1,
            max_value=20,
            value=5,
            help="Maximum number of variants to generate per compound"
        )
    
    # Submit button
    submitted = st.form_submit_button("Save Configuration", type="primary")

# Handle form submission
if submitted:
    # Validate required fields
    if not pdb_file or not resi_list:
        st.error("Please fill in all required fields marked with *")
    else:
        # Check if output directory already exists
        project_root = Path(__file__).parent.parent 
        output_path = project_root / "outputs" / output_dir_name
        if output_path.exists():
            st.error(f"Output directory '{output_dir_name}' already exists. Please choose a different name.")
        else:
            # Create temporary directory for file uploads if needed
            temp_dir = Path(tempfile.mkdtemp())
            
            # Save uploaded files
            config = {}
            
            # Always use default checkpoint path from pipeline
            config["checkpoint"] = "src/DiffSBDD/checkpoints/crossdocked_fullatom_cond.ckpt"
            
            if pdb_file:
                pdb_path = temp_dir / pdb_file.name
                with open(pdb_path, "wb") as f:
                    f.write(pdb_file.getbuffer())
                config["pdbfile"] = str(pdb_path)
            
            # Handle receptor file if uploaded
            if receptor_file:
                receptor_path = temp_dir / receptor_file.name
                with open(receptor_path, "wb") as f:
                    f.write(receptor_file.getbuffer())
                config["receptor"] = str(receptor_path)
                st.success(f"Receptor file '{receptor_file.name}' saved successfully!")
            
            # Update configuration in session state
            config.update({
                "resi_list": resi_list,
                "n_samples": n_samples,
                "sanitize": sanitize,
                "program_choice": program_choice,
                "scoring_function": scoring_function,
                "center": [center_x, center_y, center_z],
                "box_size": [size_x, size_y, size_z],
                "exhaustiveness": exhaustiveness,
                "is_selfies": is_selfies,
                "is_peptide": is_peptide,
                "top_n": top_n,
                "max_variants": max_variants,
                "num_rounds": num_rounds,
                "out_dir": str(output_path)
            })
            
            st.session_state.pipeline_config = config
            st.success("Configuration saved successfully! Proceed to the Execution page to run the pipeline.")

# Display current configuration
if st.session_state.pipeline_config is not None:
    st.header("Current Configuration")
    
    # Add stop button if pipeline is running
    if st.session_state.get("pipeline_status", {}).get("running", False):
        if st.button("Stop Pipeline", type="secondary", key="stop_button"):
            st.session_state.stop_pipeline = True
            st.warning("Stopping pipeline execution... Please wait.")
    
    st.json(st.session_state.pipeline_config) 