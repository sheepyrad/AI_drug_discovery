import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import base64
import io
import os
import json
import time
import streamlit.components.v1 as components

# Function to render molecule
def render_mol(smiles, width=400, height=300):
    """Render molecule using RDKit"""
    try:
        if smiles is None or pd.isna(smiles) or smiles == "":
            st.warning("No valid SMILES string provided")
            return None
            
        # Clean the SMILES string
        smiles = str(smiles).strip()
        
        # Try to generate the molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            st.warning(f"Could not parse SMILES: {smiles}")
            return None
            
        # Generate 2D coordinates if they don't exist
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
            
        # Create the image
        img = Draw.MolToImage(mol, size=(width, height))
        return img
    except Exception as e:
        st.error(f"Error rendering molecule: {str(e)}")
        return None

# Function to load results
def load_results(output_dir):
    """Load results from the output directory"""
    output_dir = Path(output_dir)
    results = {
        "tracking_report": None
    }
    
    # Load tracking report
    tracking_file = output_dir / "master_tracking" / "master_compound_tracking_report.csv"
    
    if tracking_file.exists():
        try:
            results["tracking_report"] = pd.read_csv(tracking_file)
            st.success("Successfully loaded tracking report.")
        except Exception as e:
            st.error(f"Error reading tracking report: {str(e)}")
            import traceback
            st.code(traceback.format_exc())
    else:
        st.warning(f"Tracking report not found at: {tracking_file}")
    
    return results

# Function to create downloadable link
def get_download_link(df, filename, text):
    """Create a download link for a dataframe"""
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">{text}</a>'
    return href

# Page configuration
st.set_page_config(
    page_title="Visualize Results",
    page_icon="🔍",
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

st.title("🔍 Visualize Pipeline Results")
st.markdown("""
    View the results from an existing output directory. Upload or select a directory to visualize results.
""")

# Initialize session state variables
if "output_dir" not in st.session_state:
    st.session_state.output_dir = None

if "results_data" not in st.session_state:
    st.session_state.results_data = None

if "selected_view" not in st.session_state:
    st.session_state.selected_view = "Summary"

# Directory selection
st.header("Select Output Directory")

# Option 1: Enter path directly
dir_path = st.text_input("Enter the path to the output directory:", placeholder="/path/to/output/directory")

# Option 2: Use directory browser
if dir_path:
    output_dir_path = Path(dir_path)
    if output_dir_path.exists() and output_dir_path.is_dir():
        st.success(f"Found directory: {output_dir_path}")
        st.session_state.output_dir = output_dir_path
    else:
        st.error(f"Directory not found: {output_dir_path}")
        st.session_state.output_dir = None

# Load results if directory is available
if st.session_state.output_dir:
    if st.button("Load Results"):
        with st.spinner("Loading results..."):
            try:
                results = load_results(st.session_state.output_dir)
                if results is None:
                    st.error("Failed to load results")
                else:
                    st.session_state.results_data = results
                    st.success("Successfully loaded results!")
            except Exception as e:
                st.error(f"Error processing results: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

# Navigation
if st.session_state.results_data:
    # Sidebar for navigation
    with st.sidebar:
        st.title("Dashboard Navigation")
        st.session_state.selected_view = st.radio(
            "Select View",
            ["Summary", "Compounds", "Variants", "Docking Results"]
        )
        
        st.divider()
        
        # Sidebar filtering options (global)
        st.subheader("Global Filters")
        df = st.session_state.results_data["tracking_report"]
        
        # Round filter (applies to all views)
        sidebar_rounds = st.multiselect(
            "Filter by Round",
            options=sorted(df["round"].unique()),
            default=sorted(df["round"].unique()),
            key="sidebar_rounds"
        )
        
        # Status filter
        sidebar_status = st.multiselect(
            "Filter by Status",
            options=sorted(df["status"].unique()),
            default=sorted(df["status"].unique()),
            key="sidebar_status"
        )
        
        # Apply global filters
        if sidebar_rounds and sidebar_status:
            filtered_df = df[df["round"].isin(sidebar_rounds) & df["status"].isin(sidebar_status)]
        else:
            filtered_df = df
            
    # Main content based on selected view
    if st.session_state.selected_view == "Summary":
        st.header("Pipeline Summary")
        
        # Summary metrics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total Compounds", len(df[df["status"] == "GENERATED"]))
        with col2:
            st.metric("Total Variants", len(df[df["status"] == "SYNTHETIZED"]))
        with col3:
            st.metric("Filtered Variants", len(df[df["status"] == "PASSFILTER"]))
        with col4:
            st.metric("Docked Compounds", len(df[df["status"] == "DOCKED"]))
        
        # Docking score distribution if available
        if "docking_score" in df.columns and df["docking_score"].notna().any():
            st.subheader("Docking Score Distribution")
            fig = px.histogram(
                df[df["docking_score"].notna()],
                x="docking_score",
                nbins=20,
                title="Distribution of Docking Scores",
                color_discrete_sequence=["#4287f5"]
            )
            fig.update_layout(bargap=0.1)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No docking scores available in the tracking report.")
    
    elif st.session_state.selected_view == "Compounds":
        st.header("Generated Compounds")
        
        # Enhanced filter controls
        filter_col1, filter_col2, filter_col3 = st.columns([2, 1, 1])
        with filter_col1:
            selected_rounds = st.multiselect(
                "Filter by Round",
                options=sorted(df["round"].unique()),
                default=sorted(df["round"].unique()),
                key="compounds_rounds"
            )
        with filter_col2:
            sort_by = st.selectbox(
                "Sort by",
                options=["compound_id", "round", "generation"],
                index=0,
                key="compounds_sort"
            )
        with filter_col3:
            sort_order = st.radio(
                "Order",
                options=["Ascending", "Descending"],
                horizontal=True,
                key="compounds_order"
            )
        
        # Filter and sort compounds
        ascending = sort_order == "Ascending"
        compounds_df = df[
            (df["status"] == "GENERATED") &
            (df["round"].isin(selected_rounds))
        ].sort_values(sort_by, ascending=ascending)
        
        # Display count
        st.info(f"Displaying {len(compounds_df)} compounds")
        
        # First show the dataframe
        st.dataframe(
            compounds_df[["compound_id", "barcode", "round", "generation", "smiles"]],
            use_container_width=True,
            hide_index=True
        )
        
        # Add download button for this filtered view
        if not compounds_df.empty:
            st.download_button(
                "Download Filtered Compounds",
                data=compounds_df.to_csv(index=False).encode('utf-8'),
                file_name="filtered_compounds.csv",
                mime="text/csv"
            )
        
        # Then show expandable elements with molecule renderings
        st.subheader("Compound Structures")
        
        # Add pagination for large datasets
        if len(compounds_df) > 10:
            compounds_per_page = st.slider("Compounds per page", 5, 20, 10, key="compounds_per_page")
            page_number = st.number_input("Page", min_value=1, max_value=max(1, len(compounds_df) // compounds_per_page + 1), step=1, key="compounds_page")
            start_idx = (page_number - 1) * compounds_per_page
            end_idx = min(start_idx + compounds_per_page, len(compounds_df))
            paginated_df = compounds_df.iloc[start_idx:end_idx]
        else:
            paginated_df = compounds_df
        
        for _, compound in paginated_df.iterrows():
            with st.expander(f"Compound {compound.get('compound_id', 'Unknown')}"):
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    # Display 2D structure
                    mol_img = render_mol(compound["smiles"])
                    if mol_img:
                        st.image(mol_img, caption="2D Structure", use_container_width=True)
                    else:
                        st.info("Could not render molecule structure")
                
                with col2:
                    # Compound details
                    details_tab, variants_tab = st.tabs(["Details", "Related Variants"])
                    
                    with details_tab:
                        st.json({
                            "Compound ID": compound.get("compound_id", ""),
                            "Barcode": compound.get("barcode", ""),
                            "Generation": compound.get("generation", ""),
                            "Round": compound.get("round", ""),
                            "SMILES": compound.get("smiles", ""),
                            "Source": compound.get("source", ""),
                            "Timestamp": compound.get("timestamp", "")
                        })
                    
                    with variants_tab:
                        # Find related variants
                        if "parent_id" in df.columns:
                            compound_id = compound.get("compound_id", "")
                            related_variants = df[df["parent_id"] == compound_id]
                            
                            if not related_variants.empty:
                                st.dataframe(
                                    related_variants[["variant_id", "status", "score"]],
                                    use_container_width=True,
                                    hide_index=True
                                )
                            else:
                                st.info("No related variants found")
                        else:
                            st.info("Parent-variant relationship not available")
    
    elif st.session_state.selected_view == "Variants":
        st.header("Synthesized Variants")
        
        # Enhanced filter controls
        filter_col1, filter_col2, filter_col3 = st.columns(3)
        with filter_col1:
            selected_status = st.multiselect(
                "Filter by Status",
                options=sorted(df["status"].unique()),
                default=[s for s in df["status"].unique() if s in ["SYNTHETIZED", "PASSFILTER"]],
                key="variants_status"
            )
        with filter_col2:
            selected_rounds_var = st.multiselect(
                "Filter by Round",
                options=sorted(df["round"].unique()),
                default=sorted(df["round"].unique()),
                key="variants_rounds"
            )
        with filter_col3:
            if "score" in df.columns:
                min_score, max_score = float(df["score"].min() if not pd.isna(df["score"]).all() else 0), float(df["score"].max() if not pd.isna(df["score"]).all() else 1)
                score_range = st.slider(
                    "Score Range",
                    min_value=min_score,
                    max_value=max_score,
                    value=(min_score, max_score),
                    key="variants_score"
                )
                score_filter = (df["score"] >= score_range[0]) & (df["score"] <= score_range[1])
            else:
                score_filter = pd.Series(True, index=df.index)
        
        # Add sorting options
        sort_col1, sort_col2 = st.columns(2)
        with sort_col1:
            sort_by_var = st.selectbox(
                "Sort by",
                options=["variant_id", "round", "generation", "score"] if "score" in df.columns else ["variant_id", "round", "generation"],
                index=0,
                key="variants_sort"
            )
        with sort_col2:
            sort_order_var = st.radio(
                "Order",
                options=["Ascending", "Descending"],
                horizontal=True,
                key="variants_order"
            )
        
        # Filter and display variants
        ascending_var = sort_order_var == "Ascending"
        variants_df = df[
            (df["status"].isin(selected_status)) &
            (df["round"].isin(selected_rounds_var)) &
            score_filter
        ].sort_values(sort_by_var, ascending=ascending_var)
        
        # Display count
        st.info(f"Displaying {len(variants_df)} variants")
        
        # Show the dataframe
        if "score" in df.columns:
            st.dataframe(
                variants_df[["variant_id", "parent_id", "round", "generation", "status", "score", "smiles"]],
                use_container_width=True,
                hide_index=True
            )
        else:
            st.dataframe(
                variants_df[["variant_id", "parent_id", "round", "generation", "status", "smiles"]],
                use_container_width=True,
                hide_index=True
            )
        
        # Add download button for this filtered view
        if not variants_df.empty:
            st.download_button(
                "Download Filtered Variants",
                data=variants_df.to_csv(index=False).encode('utf-8'),
                file_name="filtered_variants.csv",
                mime="text/csv"
            )
        
        # Variant Structures
        st.subheader("Variant Structures")
        
        # Add pagination for large datasets
        if len(variants_df) > 10:
            variants_per_page = st.slider("Variants per page", 5, 20, 10, key="variants_per_page")
            page_number = st.number_input("Page", min_value=1, max_value=max(1, len(variants_df) // variants_per_page + 1), step=1, key="variants_page")
            start_idx = (page_number - 1) * variants_per_page
            end_idx = min(start_idx + variants_per_page, len(variants_df))
            paginated_var_df = variants_df.iloc[start_idx:end_idx]
        else:
            paginated_var_df = variants_df
        
        for _, variant in paginated_var_df.iterrows():
            with st.expander(f"Variant {variant.get('variant_id', 'Unknown')}"):
                col1, col2 = st.columns([1, 2])
                
                with col1:
                    # Display 2D structure
                    mol_img = render_mol(variant["smiles"])
                    if mol_img:
                        st.image(mol_img, caption="2D Structure", use_container_width=True)
                    else:
                        st.info("Could not render molecule structure")
                    
                    # If parent molecule exists, display it for comparison
                    if "parent_id" in variant and "source_smiles" in variant:
                        st.subheader("Parent Structure")
                        parent_mol_img = render_mol(variant["source_smiles"])
                        if parent_mol_img:
                            st.image(parent_mol_img, caption="Parent Structure", use_container_width=True)
                
                with col2:
                    # Variant details
                    details = {
                        "Variant ID": variant.get("variant_id", ""),
                        "Parent Compound": variant.get("parent_id", ""),
                        "Status": variant.get("status", ""),
                        "Generation": variant.get("generation", ""),
                        "Round": variant.get("round", ""),
                        "Source": variant.get("source", ""),
                        "Timestamp": variant.get("timestamp", "")
                    }
                    if "score" in variant:
                        details["Retrosynthesis Score"] = variant["score"]
                    st.json(details)
                    
                    # Link to parent
                    parent_id = variant.get("parent_id", "")
                    if parent_id and parent_id in df["compound_id"].values:
                        parent_smiles = df[df["compound_id"] == parent_id]["smiles"].values[0]
                        st.markdown(f"**Parent SMILES**: `{parent_smiles}`")
    
    elif st.session_state.selected_view == "Docking Results":
        st.header("Docking Results")
        
        if "docking_score" not in df.columns:
            st.info("No docking scores available in the tracking report.")
        else:
            # Enhanced filter controls for docking
            filter_col1, filter_col2 = st.columns(2)
            with filter_col1:
                selected_rounds_dock = st.multiselect(
                    "Filter by Round",
                    options=sorted(df["round"].unique()),
                    default=sorted(df["round"].unique()),
                    key="docking_rounds"
                )
            with filter_col2:
                min_score, max_score = float(df["docking_score"].min() if not pd.isna(df["docking_score"]).all() else 0), float(df["docking_score"].max() if not pd.isna(df["docking_score"]).all() else 0)
                score_range_dock = st.slider(
                    "Docking Score Range",
                    min_value=min_score,
                    max_value=max_score,
                    value=(min_score, max_score),
                    key="docking_score_range"
                )
            
            # Filter docked compounds
            docked_df = df[
                (df["status"] == "DOCKED") &
                (df["round"].isin(selected_rounds_dock)) &
                (df["docking_score"] >= score_range_dock[0]) &
                (df["docking_score"] <= score_range_dock[1])
            ].sort_values("docking_score")
            
            if docked_df.empty:
                st.info("No docked compounds match the current filter criteria.")
            else:
                # Docking statistics
                stats_cols = st.columns(4)
                with stats_cols[0]:
                    st.metric("Best Score", f"{docked_df['docking_score'].min():.2f}")
                with stats_cols[1]:
                    st.metric("Average Score", f"{docked_df['docking_score'].mean():.2f}")
                with stats_cols[2]:
                    st.metric("Median Score", f"{docked_df['docking_score'].median():.2f}")
                with stats_cols[3]:
                    st.metric("Total Docked", len(docked_df))
                
                # Score distribution
                st.subheader("Score Distribution")
                fig = px.scatter(
                    docked_df,
                    x="round",
                    y="docking_score",
                    color="docking_score",
                    hover_data=["compound_id", "smiles"],
                    title="Docking Scores by Round",
                    color_continuous_scale="viridis"
                )
                st.plotly_chart(fig, use_container_width=True)
                
                # Full docking results table
                st.subheader("All Docking Results")
                st.dataframe(
                    docked_df[["compound_id", "variant_id", "round", "docking_score", "status", "smiles"]],
                    use_container_width=True,
                    hide_index=True
                )
                
                # Add download button for this filtered view
                st.download_button(
                    "Download Filtered Docking Results",
                    data=docked_df.to_csv(index=False).encode('utf-8'),
                    file_name="filtered_docking_results.csv",
                    mime="text/csv"
                )
                
                # Top results
                st.subheader("Top Docking Results")
                
                # Add pagination for large datasets
                if len(docked_df) > 10:
                    top_n = st.slider("Show top N results", 5, min(30, len(docked_df)), 10, key="docking_top_n")
                    top_results = docked_df.head(top_n)
                else:
                    top_results = docked_df
                
                for _, result in top_results.iterrows():
                    # Use variant_id if available, otherwise fallback to compound_id
                    display_id = result.get('variant_id', result.get('compound_id', 'Unknown'))
                    
                    with st.expander(f"{display_id} (Score: {result.get('docking_score', 0):.2f})"):
                        col1, col2 = st.columns([1, 2])
                        with col1:
                            mol_img = render_mol(result["smiles"])
                            if mol_img:
                                st.image(mol_img, caption="2D Structure", use_container_width=True)
                            else:
                                st.info("Could not render molecule structure")
                        with col2:
                            st.json({
                                "Compound ID": result.get("compound_id", ""),
                                "Variant ID": result.get("variant_id", ""),
                                "Barcode": result.get("barcode", ""),
                                "Round": result.get("round", ""),
                                "Docking Score": result.get("docking_score", ""),
                                "Best Pose": result.get("best_pose", ""),
                                "SMILES": result.get("smiles", "")
                            })
                            
                            # Download options for 3D files
                            if "barcode" in result and "best_pose" in result:
                                # Path to the docked ligand pose file
                                variant_dir = Path(st.session_state.output_dir) / f"round_{result['round']}" / "docking_results" / f"variant_{result['barcode']}"
                                
                                # Parse best pose value (handle decimal values by converting to integer)
                                try:
                                    if pd.notna(result["best_pose"]):
                                        best_pose_num = int(float(result["best_pose"]))
                                        pose_file = variant_dir / f"pose_{best_pose_num}.pdbqt"
                                    else:
                                        pose_file = None
                                except (ValueError, TypeError):
                                    # If conversion fails, look for any pose files
                                    pose_file = None
                                
                                # Get receptor file path from the input directory
                                receptor_file = Path(st.session_state.output_dir).parent / "input" / "NS5_test.pdbqt"
                                
                                # Check if any pose file exists
                                if pose_file is not None and pose_file.exists():
                                    # Download buttons for the best pose and receptor
                                    st.download_button(
                                        "Download Best Pose",
                                        data=open(pose_file, 'rb').read(),
                                        file_name=pose_file.name,
                                        mime="application/octet-stream",
                                        key=f"dl_best_pose_{display_id}"
                                    )
                                    
                                    if receptor_file.exists():
                                        st.download_button(
                                            "Download Receptor File",
                                            data=open(receptor_file, 'rb').read(),
                                            file_name=receptor_file.name,
                                            mime="application/octet-stream",
                                            key=f"dl_receptor_best_{display_id}"
                                        )
                                else:
                                    # Always try to find all available pose files in the variant directory
                                    if variant_dir.exists():
                                        pose_files = list(variant_dir.glob("pose_*.pdbqt"))
                                        if pose_files:
                                            st.info(f"Found {len(pose_files)} pose files")
                                            
                                            # Extract pose numbers and sort them
                                            pose_numbers = []
                                            for p in pose_files:
                                                try:
                                                    # Extract number from pose_X.pdbqt
                                                    pose_num = int(p.stem.split('_')[1])
                                                    pose_numbers.append((pose_num, p))
                                                except (ValueError, IndexError):
                                                    continue
                                            
                                            # Sort poses by number
                                            pose_numbers.sort()
                                            sorted_poses = [p[1] for p in pose_numbers]
                                            
                                            # Create a selection for the poses
                                            pose_options = [f"Pose {p[0]}" for p in pose_numbers]
                                            selected_pose_idx = st.selectbox(
                                                "Select pose", 
                                                range(len(pose_options)),
                                                format_func=lambda i: pose_options[i],
                                                key=f"pose_select_{display_id}"
                                            )
                                            
                                            selected_file = sorted_poses[selected_pose_idx]
                                            
                                            # Download buttons for the selected pose and receptor
                                            st.download_button(
                                                "Download Selected Pose",
                                                data=open(selected_file, 'rb').read(),
                                                file_name=selected_file.name,
                                                mime="application/octet-stream",
                                                key=f"dl_pose_{display_id}"
                                            )
                                            
                                            if receptor_file.exists():
                                                st.download_button(
                                                    "Download Receptor File",
                                                    data=open(receptor_file, 'rb').read(),
                                                    file_name=receptor_file.name,
                                                    mime="application/octet-stream",
                                                    key=f"dl_receptor_{display_id}"
                                                )
                                        else:
                                            st.info(f"No pose files found in {variant_dir}")
                                    else:
                                        st.info(f"Variant directory not found: {variant_dir}")
                            else:
                                st.info("Barcode or best pose information missing for 3D structure visualization")
    
    # Export options
    st.divider()
    st.subheader("Export Options")
    
    # Add download button for complete dataset
    if st.button("Export All Results"):
        st.download_button(
            "📥 Download Complete Dataset",
            data=df.to_csv(index=False).encode('utf-8'),
            file_name="all_results.csv",
            mime="text/csv"
        )
else:
    st.info("Please select an output directory and load results to view visualizations.") 