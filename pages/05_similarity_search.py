import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import Draw
import plotly.express as px
import plotly.graph_objects as go
import os
from PIL import Image
import io

# --- Configuration ---
PAGE_TITLE = "Similarity Search"

st.set_page_config(page_title=PAGE_TITLE, layout="wide")
st.title(PAGE_TITLE)
st.markdown("""
Perform Tanimoto similarity searches against compounds in the master tracking report.
Upload a CSV file or provide the path to a local CSV file.
Enter a SMILES string for the query molecule and select filters to analyze the results.
""")

# --- Data Loading ---
@st.cache_data
def load_data(data_source):
    """Loads the master compound tracking report from path or uploaded file."""
    source_name = ""
    try:
        if isinstance(data_source, str): # Check if it's a path
            if not os.path.exists(data_source):
                st.sidebar.error(f"Error: File not found at path: {data_source}") # Show error in sidebar
                return None, data_source # Return path even on error
            source_name = data_source
            df = pd.read_csv(data_source)
        elif data_source is not None: # Assume it's an uploaded file object
            source_name = data_source.name
            df = pd.read_csv(data_source)
        else:
            return None, "No data source provided"

        # Basic data cleaning/validation
        required_cols = ['compound_id', 'round', 'status', 'smiles']
        if not all(col in df.columns for col in required_cols):
            # Use sidebar for persistent errors related to loading
            st.sidebar.error(f"CSV ({source_name}) missing required columns: {required_cols}. Found: {list(df.columns)}")
            return None, source_name
        df = df.dropna(subset=['smiles']) # Drop rows with missing SMILES
        df['round'] = pd.to_numeric(df['round'], errors='coerce').astype('Int64') # Ensure round is integer
        # Use sidebar for success message
        st.sidebar.success(f"Loaded {len(df)} compounds from `{source_name}`")
        return df, source_name
    except Exception as e:
        st.sidebar.error(f"Error loading data from {source_name or data_source}: {e}")
        return None, source_name # Return source_name even on error for display

# --- RDKit Functions ---
def calculate_tanimoto_from_smiles(query_smiles, target_df):
    """Calculate Tanimoto similarity between query SMILES and all compounds in the dataframe."""
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        return None, "Invalid SMILES string for query molecule."

    query_fp = Chem.RDKFingerprint(query_mol)
    if query_fp is None:
        return None, "Could not generate fingerprint for query molecule."

    similarities = []
    skipped_count = 0
    for idx, row in target_df.iterrows():
        try:
            if pd.isna(row['smiles']):
                skipped_count += 1
                continue

            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is None:
                similarities.append({
                    'compound_id': row.get('compound_id', f'Index_{idx}'),
                    'round': row.get('round', None),
                    'status': row.get('status', 'N/A'),
                    'tanimoto_similarity': np.nan, # Use NaN for invalid SMILES
                    'smiles': row['smiles']
                })
                skipped_count += 1
                continue

            fp = Chem.RDKFingerprint(mol)
            if fp is None:
                 similarities.append({
                    'compound_id': row.get('compound_id', f'Index_{idx}'),
                    'round': row.get('round', None),
                    'status': row.get('status', 'N/A'),
                    'tanimoto_similarity': np.nan, # Use NaN if fingerprint fails
                    'smiles': row['smiles']
                })
                 skipped_count += 1
                 continue

            sim = DataStructs.TanimotoSimilarity(query_fp, fp)

            similarities.append({
                'compound_id': row.get('compound_id', f'Index_{idx}'),
                'round': row.get('round', None),
                'status': row.get('status', 'N/A'),
                'tanimoto_similarity': sim,
                'smiles': row['smiles']
            })
        except Exception as e:
            # Use st.warning for non-critical processing errors during calculation
            st.warning(f"Error processing compound {row.get('compound_id', f'Index_{idx}')}: {e}")
            similarities.append({
                'compound_id': row.get('compound_id', f'Index_{idx}'),
                'round': row.get('round', None),
                'status': row.get('status', 'N/A'),
                'tanimoto_similarity': np.nan, # Use NaN on error
                'smiles': row.get('smiles', 'Error')
            })
            skipped_count += 1

    if skipped_count > 0:
        st.warning(f"Skipped {skipped_count} compounds during similarity calculation due to missing/invalid SMILES or processing errors.")

    return pd.DataFrame(similarities), None

def draw_molecule(smiles):
    """Draw a molecule from SMILES string and return PIL image"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        img = Draw.MolToImage(mol, size=(300, 300))
        return img
    except Exception as e:
        st.error(f"Error generating molecule image: {e}")
        return None

# --- Plotting Functions ---
def plot_similarity_by_round(similarity_df, title="Tanimoto Similarity by Round"):
    """Plot Tanimoto similarity grouped by round"""
    filtered_df = similarity_df.dropna(subset=['tanimoto_similarity']).copy()
    if filtered_df.empty:
        st.warning("No valid similarity data to plot by round.")
        return go.Figure()
    # Ensure 'round' exists and handle potential missing values before conversion
    if 'round' not in filtered_df.columns:
        st.warning("'round' column missing for round plot.")
        return go.Figure()
    filtered_df = filtered_df.dropna(subset=['round'])
    filtered_df['round'] = filtered_df['round'].astype(str) # Plotly box plot needs string category
    if filtered_df.empty:
        st.warning("No valid round data after filtering NAs.")
        return go.Figure()

    fig = px.box(filtered_df, x='round', y='tanimoto_similarity',
                 title=title,
                 labels={'tanimoto_similarity': 'Tanimoto Similarity', 'round': 'Round'},
                 color='round',
                 hover_data=['compound_id', 'status', 'smiles'])

    fig.add_trace(go.Scatter(x=filtered_df['round'], y=filtered_df['tanimoto_similarity'],
                             mode='markers', name='Compounds',
                             marker=dict(opacity=0.5, size=5),
                             customdata=filtered_df[['compound_id', 'status', 'smiles']],
                             hovertemplate='<b>Round</b>: %{x}<br>' +
                                           '<b>Similarity</b>: %{y:.3f}<br>' +
                                           '<b>ID</b>: %{customdata[0]}<br>' +
                                           '<b>Status</b>: %{customdata[1]}<br>' +
                                           '<b>SMILES</b>: %{customdata[2]}<extra></extra>',
                             showlegend=False))

    fig.update_layout(height=600)
    return fig

def plot_similarity_by_status(similarity_df, title="Tanimoto Similarity by Status"):
    """Plot Tanimoto similarity grouped by status"""
    filtered_df = similarity_df.dropna(subset=['tanimoto_similarity']).copy()
    if filtered_df.empty:
        st.warning("No valid similarity data to plot by status.")
        return go.Figure()
    # Ensure 'status' exists
    if 'status' not in filtered_df.columns:
        st.warning("'status' column missing for status plot.")
        return go.Figure()
    filtered_df = filtered_df.dropna(subset=['status'])
    if filtered_df.empty:
        st.warning("No valid status data after filtering NAs.")
        return go.Figure()

    fig = px.box(filtered_df, x='status', y='tanimoto_similarity',
                 title=title,
                 labels={'tanimoto_similarity': 'Tanimoto Similarity', 'status': 'Status'},
                 color='status',
                 hover_data=['compound_id', 'round', 'smiles'])

    fig.add_trace(go.Scatter(x=filtered_df['status'], y=filtered_df['tanimoto_similarity'],
                             mode='markers', name='Compounds',
                             marker=dict(opacity=0.5, size=5),
                             customdata=filtered_df[['compound_id', 'round', 'smiles']],
                             hovertemplate='<b>Status</b>: %{x}<br>' +
                                           '<b>Similarity</b>: %{y:.3f}<br>' +
                                           '<b>ID</b>: %{customdata[0]}<br>' +
                                           '<b>Round</b>: %{customdata[1]}<br>' +
                                           '<b>SMILES</b>: %{customdata[2]}<extra></extra>',
                             showlegend=False))

    fig.update_layout(height=600)
    return fig

def plot_similarity_heatmap(similarity_df, title="Average Tanimoto Similarity Heatmap"):
    """Plot Tanimoto similarity as a heatmap by round and status"""
    filtered_df = similarity_df.dropna(subset=['tanimoto_similarity']).copy()
    if filtered_df.empty or 'round' not in filtered_df.columns or 'status' not in filtered_df.columns:
        st.warning("Insufficient data for heatmap (requires valid similarity, round, and status).")
        return go.Figure()

    # Ensure round is numeric for pivot table, handle potential NaNs
    filtered_df['round'] = pd.to_numeric(filtered_df['round'], errors='coerce')
    filtered_df = filtered_df.dropna(subset=['round', 'status']) # Need both for pivot
    filtered_df['round'] = filtered_df['round'].astype(int)

    if filtered_df.empty:
        st.warning("No valid numeric round/status data for heatmap after filtering NAs.")
        return go.Figure()

    try:
        pivot_df = filtered_df.pivot_table(
            values='tanimoto_similarity',
            index='round',
            columns='status',
            aggfunc='mean'
        )
        if pivot_df.empty:
            st.warning("Pivot table is empty. Cannot generate heatmap.")
            return go.Figure()

        fig = px.imshow(pivot_df,
                        labels=dict(x="Status", y="Round", color="Avg. Tanimoto Similarity"),
                        title=title,
                        color_continuous_scale="Viridis",
                        text_auto=".2f") # Display values on heatmap

        fig.update_layout(height=600)
        return fig
    except Exception as e:
        st.error(f"Error creating heatmap pivot table: {e}")
        return go.Figure()

def plot_top_similar_compounds(similarity_df, n=10, title="Top N Similar Compounds"):
    """Plot the top N most similar compounds"""
    filtered_df = similarity_df.dropna(subset=['tanimoto_similarity']).copy()
    if filtered_df.empty:
        st.warning("No valid similarity data to find top compounds.")
        return go.Figure(), pd.DataFrame()

    top_df = filtered_df.sort_values('tanimoto_similarity', ascending=False).head(n)
    if top_df.empty:
        st.warning(f"No compounds found after sorting for top {n}.")
        return go.Figure(), pd.DataFrame()

    top_df['similarity_rank'] = range(1, len(top_df) + 1)

    # Ensure required columns exist for hover data
    hover_cols = ['round', 'smiles', 'similarity_rank']
    available_hover_cols = [col for col in hover_cols if col in top_df.columns]

    fig = px.bar(top_df, x='compound_id', y='tanimoto_similarity',
                 color='status' if 'status' in top_df.columns else None, # Handle missing status
                 hover_data=available_hover_cols,
                 title=f"{title} (Top {min(n, len(top_df))})", # Adjust title based on actual count
                 labels={'tanimoto_similarity': 'Tanimoto Similarity', 'compound_id': 'Compound ID'})

    fig.update_layout(height=600, xaxis_tickangle=-45)
    # Define cols needed for table output
    table_cols = ['similarity_rank', 'compound_id', 'tanimoto_similarity', 'smiles', 'round', 'status']
    available_table_cols = [col for col in table_cols if col in top_df.columns]
    return fig, top_df[available_table_cols]

# --- UI Elements --- #

# Initialize session state variables if they don't exist
if 'data_source_name' not in st.session_state:
    st.session_state.data_source_name = "No data loaded"
if 'df' not in st.session_state:
    st.session_state.df = None
if 'uploaded_file_state' not in st.session_state:
    st.session_state.uploaded_file_state = None
if 'file_path_state' not in st.session_state:
    st.session_state.file_path_state = ""

# Configuration Expander in Main Area
with st.expander("Configuration & Query", expanded=True):
    st.subheader("1. Data Source")
    uploaded_file = st.file_uploader("Upload Master Tracking Report (CSV)", type="csv", key="file_uploader")
    file_path_input = st.text_input("Or Enter Path to CSV File", placeholder="e.g., outputs/NS5/master_tracking/report.csv", key="path_input")

    # Check if data source has changed or needs loading
    load_trigger = False
    current_data_input = uploaded_file if uploaded_file is not None else file_path_input
    previous_data_input_state = st.session_state.uploaded_file_state if st.session_state.uploaded_file_state is not None else st.session_state.file_path_state

    # Trigger loading if uploaded file changes, or if path input changes and is not empty
    if uploaded_file is not None and uploaded_file != st.session_state.uploaded_file_state:
        st.session_state.file_path_state = "" # Clear path if upload is used
        file_path_input = "" # Reset widget state visually if needed
        st.session_state.uploaded_file_state = uploaded_file
        load_trigger = True
    elif file_path_input != "" and file_path_input != st.session_state.file_path_state:
        st.session_state.uploaded_file_state = None # Clear upload if path is used
        # uploaded_file = None # Reset widget state maybe? Might clear user selection unexpectedly.
        st.session_state.file_path_state = file_path_input
        load_trigger = True
    elif uploaded_file is None and st.session_state.uploaded_file_state is not None: # File removed
        st.session_state.uploaded_file_state = None
        load_trigger = True # Trigger reload (will likely result in df=None)
    elif file_path_input == "" and st.session_state.file_path_state != "": # Path cleared
        st.session_state.file_path_state = ""
        load_trigger = True # Trigger reload (will likely result in df=None)


    if load_trigger:
        data_source_to_load = st.session_state.uploaded_file_state if st.session_state.uploaded_file_state is not None else st.session_state.file_path_state
        if data_source_to_load: # Check if there is actually something to load
            st.session_state.df, st.session_state.data_source_name = load_data(data_source_to_load)
        else: # Both are empty/None
             st.session_state.df, st.session_state.data_source_name = None, "No data loaded"
        st.rerun() # Rerun to update UI based on new df state

    # Use loaded data from session state
    df = st.session_state.df
    data_source_name = st.session_state.data_source_name

    # Display data source info in the sidebar regardless
    st.sidebar.markdown("---")
    st.sidebar.info(f"Current Data Source: `{data_source_name}`")

    if df is not None:
        st.subheader("2. Query and Filters")
        query_smiles = st.text_input("Query SMILES", placeholder="Enter SMILES string (e.g., COc1ccccc1)", key="smiles_input")

        col1, col2 = st.columns(2)
        with col1:
            # Ensure status column exists before creating widget
            if 'status' in df.columns:
                 available_statuses = sorted(df['status'].astype(str).unique()) # Handle potential non-string types
                 selected_statuses = st.multiselect("Filter by Status", options=available_statuses, default=available_statuses, key="status_filter")
            else:
                st.warning("'status' column not found in data. Cannot filter by status.")
                selected_statuses = []
        with col2:
            # Ensure round column exists and handle potential errors
            if 'round' in df.columns:
                try:
                    available_rounds = sorted(df['round'].dropna().unique().astype(int))
                    selected_rounds = st.multiselect("Filter by Round", options=available_rounds, default=available_rounds, key="round_filter")
                except Exception as e:
                     st.warning(f"Could not process 'round' column for filtering: {e}. Ensure it contains numeric values.")
                     selected_rounds = []
            else:
                st.warning("'round' column not found in data. Cannot filter by round.")
                selected_rounds = []

        st.subheader("3. Analysis Options")
        plot_type = st.selectbox(
            "Plot Type",
            options=['By Round', 'By Status', 'Heatmap', 'Top N Similar'],
            index=0,
            key="plot_type_select"
        )

        top_n = 10
        if plot_type == 'Top N Similar':
            top_n = st.slider("Number of Top Compounds (N)", min_value=5, max_value=50, value=10, step=5, key="top_n_slider")

        calculate_button = st.button("Calculate Similarity", type="primary", use_container_width=True, key="calc_button")

    else:
        # If df is None after checking inputs, show info here
        if not st.session_state.uploaded_file_state and not st.session_state.file_path_state:
             st.info("Please upload a CSV file or provide a path above to load compound data.")
        # Error messages for failed loading are now shown in the sidebar by load_data

        calculate_button = False # Disable calculation if no data
        query_smiles = "" # Reset query smiles if no data
        selected_statuses = []
        selected_rounds = []
        plot_type = 'By Round'
        top_n = 10

# --- Calculation and Display Logic (Main Area) --- #

# This block runs only if the button was clicked in the expander AND df is loaded
if calculate_button and df is not None:
    if not query_smiles:
        st.warning("Please enter a query SMILES string in the configuration section.")
        st.stop()

    st.header("Results")

    # Display Query Molecule
    st.subheader("Query Molecule")
    query_mol_img = draw_molecule(query_smiles)
    if query_mol_img:
        st.image(query_mol_img, caption=f"Query: {query_smiles}")
    else:
        st.error("Invalid Query SMILES string. Please check and try again.")
        st.stop() # Stop execution if SMILES is invalid

    # Filter DataFrame based on selections
    filtered_df_calc = df.copy()
    if selected_statuses and 'status' in filtered_df_calc.columns: # Check if list is not empty and column exists
         filtered_df_calc = filtered_df_calc[filtered_df_calc['status'].isin(selected_statuses)]
    if selected_rounds and 'round' in filtered_df_calc.columns: # Check if list is not empty and column exists
         # Ensure rounds are compared correctly (e.g., filtering ints)
         try:
             numeric_rounds = pd.to_numeric(filtered_df_calc['round'], errors='coerce')
             filtered_df_calc = filtered_df_calc[numeric_rounds.isin(selected_rounds)]
         except Exception as e:
              st.warning(f"Could not apply round filter: {e}")

    st.write(f"Comparing against {len(filtered_df_calc)} compounds based on filters.")

    if filtered_df_calc.empty:
         st.warning("No compounds match the selected filters for calculation.")
    else:
        # Calculate Similarity
        with st.spinner("Calculating Tanimoto similarities..."):
            similarity_df, error = calculate_tanimoto_from_smiles(query_smiles, filtered_df_calc)

        if error:
            st.error(f"Error during calculation: {error}")
        elif similarity_df is None or similarity_df.empty:
             st.warning("Calculation returned no results or failed for all filtered compounds.")
        else:
            valid_similarities = similarity_df.dropna(subset=['tanimoto_similarity'])

            if valid_similarities.empty:
                 st.warning("No valid similarity scores could be calculated for the filtered compounds.")
            else:
                # Display Summary Statistics
                st.subheader("Summary Statistics")
                stats_cols = st.columns(4)
                stats_cols[0].metric("Compounds Analyzed", len(valid_similarities))
                stats_cols[1].metric("Average Similarity", f"{valid_similarities['tanimoto_similarity'].mean():.4f}")
                stats_cols[2].metric("Median Similarity", f"{valid_similarities['tanimoto_similarity'].median():.4f}")
                stats_cols[3].metric("Max Similarity", f"{valid_similarities['tanimoto_similarity'].max():.4f}")

                # Display Plot
                st.subheader(f"Plot: {plot_type}")
                fig = None
                top_df_table = None

                try:
                    if plot_type == 'By Round':
                        fig = plot_similarity_by_round(similarity_df)
                    elif plot_type == 'By Status':
                        fig = plot_similarity_by_status(similarity_df)
                    elif plot_type == 'Heatmap':
                         fig = plot_similarity_heatmap(similarity_df)
                    elif plot_type == 'Top N Similar':
                        fig, top_df_table = plot_top_similar_compounds(similarity_df, n=top_n)

                    if fig is not None and len(fig.data) > 0: # Check if figure has data
                         st.plotly_chart(fig, use_container_width=True)
                    elif fig is not None: # Figure exists but is empty
                         st.info("Plot generated, but contains no data based on current filters/results or selected plot type requires missing columns.")
                    else:
                         st.warning(f"Could not generate plot for type: {plot_type}")

                except Exception as e:
                    st.error(f"An error occurred during plot generation: {e}")

                # Display Top N Table (if applicable)
                if top_df_table is not None and not top_df_table.empty:
                    st.subheader(f"Top {len(top_df_table)} Similar Compounds Table")
                    st.dataframe(top_df_table.reset_index(drop=True), use_container_width=True) # Reset index for cleaner display

                # Display Most Similar Compound Details
                st.subheader("Most Similar Compound")
                try:
                    most_similar = valid_similarities.loc[valid_similarities['tanimoto_similarity'].idxmax()]
                    cols = st.columns([1, 3])
                    with cols[0]:
                        st.metric("Similarity", f"{most_similar.get('tanimoto_similarity', 'N/A'):.4f}")
                        st.write(f"**ID:** {most_similar.get('compound_id', 'N/A')}")
                        # Handle potential missing columns gracefully
                        round_val = most_similar.get('round', 'N/A')
                        status_val = most_similar.get('status', 'N/A')
                        smiles_val = most_similar.get('smiles', 'N/A')
                        st.write(f"**Round:** {round_val}")
                        st.write(f"**Status:** {status_val}")
                        st.write(f"**SMILES:** `{smiles_val}`")
                    with cols[1]:
                        if smiles_val != 'N/A':
                            most_similar_img = draw_molecule(smiles_val)
                            if most_similar_img:
                                st.image(most_similar_img, caption=f"Most Similar: {most_similar.get('compound_id', 'N/A')}")
                            else:
                                st.warning("Could not draw molecule for the most similar compound.")
                        else:
                             st.warning("SMILES string missing for the most similar compound.")
                except Exception as e:
                     st.error(f"Error displaying most similar compound details: {e}") 