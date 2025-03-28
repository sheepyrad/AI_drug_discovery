import streamlit as st
import time
from pathlib import Path
import threading
import queue
import logging
from datetime import datetime
import os
import sys
import streamlit.components.v1 as components
import json

# Import the pipeline function
from pipeline_quick_multiround import main as pipeline_main

# Configure logging for the execution page
logger = logging.getLogger(__name__)

st.set_page_config(
    page_title="Pipeline Execution",
    page_icon="🚀",
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

# Initialize session state variables if they don't exist
if "pipeline_status" not in st.session_state:
    st.session_state.pipeline_status = {
        "running": False,
        "current_round": 0,
        "total_rounds": 0,
        "progress": 0,
        "last_log_position": 0  # Track last read position in log file
    }

if "pipeline_thread" not in st.session_state:
    st.session_state.pipeline_thread = None

if "stop_event" not in st.session_state:
    st.session_state.stop_event = threading.Event()

if "log_container" not in st.session_state:
    st.session_state.log_container = ""

if "last_refresh" not in st.session_state:
    st.session_state.last_refresh = 0

def read_log_file(log_path):
    """Read a log file with error handling and position tracking"""
    try:
        if not Path(log_path).exists():
            return None
            
        with open(log_path, 'r') as f:
            # Seek to last read position
            f.seek(st.session_state.pipeline_status["last_log_position"])
            # Read new content
            new_content = f.read()
            # Update last read position
            st.session_state.pipeline_status["last_log_position"] = f.tell()
            return new_content
    except Exception as e:
        st.error(f"Error reading log file {log_path}: {str(e)}")
        return None

def read_json_file(file_path):
    """Read and parse a JSON file with error handling"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        st.error(f"Error reading JSON file {file_path}: {str(e)}")
        return None

def update_progress_from_log(log_content):
    """Update progress based on log content"""
    try:
        # Look for round information in log
        if "STARTING ROUND" in log_content:
            current_round = int(log_content.split("STARTING ROUND")[1].split("/")[0])
            total_rounds = int(log_content.split("STARTING ROUND")[1].split("/")[1].split()[0])
            st.session_state.pipeline_status["current_round"] = current_round
            st.session_state.pipeline_status["total_rounds"] = total_rounds
            # Calculate progress as percentage of rounds completed
            progress = ((current_round - 1) / total_rounds) * 100
            st.session_state.pipeline_status["progress"] = progress
        
        # Update progress based on specific pipeline stages
        if "Running ligand generation" in log_content:
            stage_progress = 20
        elif "Running retrosynthesis" in log_content:
            stage_progress = 40
        elif "Starting batch filtering" in log_content:
            stage_progress = 60
        elif "Starting docking" in log_content:
            stage_progress = 80
        elif "Pipeline completed successfully" in log_content:
            stage_progress = 100
            st.session_state.pipeline_status["running"] = False
        
        # Update progress if we found a stage
        if 'stage_progress' in locals():
            current_round = st.session_state.pipeline_status["current_round"]
            total_rounds = st.session_state.pipeline_status["total_rounds"]
            if total_rounds > 0:
                base_progress = ((current_round - 1) / total_rounds) * 100
                round_progress = (stage_progress / 100) * (100 / total_rounds)
                st.session_state.pipeline_status["progress"] = min(base_progress + round_progress, 100)
    
    except Exception as e:
        logger.warning(f"Error updating progress: {e}")

def create_auto_scrolling_text_area(content, height=400):
    """Create an auto-scrolling text area using HTML and JavaScript with syntax highlighting"""
    # Escape HTML special characters
    content = content.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    
    # Add syntax highlighting for common log patterns
    content = content.replace("ERROR", '<span style="color: #ff6b6b;">ERROR</span>')
    content = content.replace("WARNING", '<span style="color: #ffd93d;">WARNING</span>')
    content = content.replace("INFO", '<span style="color: #6bff6b;">INFO</span>')
    content = content.replace("DEBUG", '<span style="color: #6b6bff;">DEBUG</span>')
    
    # Add syntax highlighting for pipeline stages
    stages = [
        "Running ligand generation",
        "Running retrosynthesis",
        "Starting batch filtering",
        "Starting docking",
        "Pipeline completed successfully",
        "STARTING ROUND",
        "COMPLETED ROUND"
    ]
    
    for stage in stages:
        content = content.replace(stage, f'<span style="color: #ffa500;">{stage}</span>')
    
    html = f"""
        <div style="
            height: {height}px;
            overflow-y: auto;
            border: 1px solid #2e2e2e;
            border-radius: 4px;
            padding: 12px;
            background-color: #1e1e1e;
            color: #d4d4d4;
            font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
            font-size: 14px;
            line-height: 1.5;
            margin-bottom: 20px;
        ">
            <pre id="log-content" style="margin: 0; white-space: pre-wrap; padding-bottom: 20px;">{content}</pre>
        </div>
        <script>
            // Function to scroll to bottom
            function scrollToBottom() {{
                var element = document.getElementById('log-content');
                var container = element.parentElement;
                container.scrollTop = container.scrollHeight;
                // Add a small delay to ensure the scroll happens after content is rendered
                setTimeout(() => {{
                    container.scrollTop = container.scrollHeight;
                }}, 100);
            }}
            
            // Initial scroll
            scrollToBottom();
            
            // Set up a mutation observer to watch for content changes
            var observer = new MutationObserver(function(mutations) {{
                scrollToBottom();
            }});
            
            observer.observe(document.getElementById('log-content'), {{
                childList: true,
                characterData: true,
                subtree: true
            }});
            
            // Also scroll on window resize
            window.addEventListener('resize', scrollToBottom);
        </script>
    """
    return components.html(html, height=height)

st.title("▶️ Pipeline Execution")

# Check if configuration exists
if not st.session_state.get("pipeline_config"):
    st.error("Please configure the pipeline parameters first!")
    st.stop()

def run_pipeline(config, status_dict, stop_event):
    """Function to run the pipeline"""
    try:
        # Update status using the passed dictionary
        status_dict["running"] = True
        status_dict["current_round"] = 1
        status_dict["total_rounds"] = config["num_rounds"]
        status_dict["progress"] = 0
        
        # Check for receptor file
        if "receptor" not in config or not config["receptor"]:
            logger.warning("No receptor file provided. Docking steps will be skipped or may fail.")
            # You could either set a default receptor or handle this in the pipeline
            
        # Run the pipeline with stop check
        pipeline_main(
            out_dir=config["out_dir"],
            checkpoint=config["checkpoint"],
            pdbfile=config["pdbfile"],
            resi_list=config["resi_list"],
            n_samples=config["n_samples"],
            sanitize=config["sanitize"],
            receptor=config.get("receptor"),
            program_choice=config["program_choice"],
            scoring_function=config["scoring_function"],
            center=config["center"],
            box_size=config["box_size"],
            exhaustiveness=config["exhaustiveness"],
            is_selfies=config["is_selfies"],
            is_peptide=config["is_peptide"],
            top_n=config["top_n"],
            max_variants=config["max_variants"],
            num_rounds=config["num_rounds"],
            stop_flag=status_dict  # Pass the status dict for stop checking
        )
        
        status_dict["running"] = False
        status_dict["progress"] = 100
        
    except Exception as e:
        logger.error(f"Pipeline execution failed: {str(e)}")
        status_dict["running"] = False
        status_dict["error"] = str(e)
    finally:
        # Clear the stop event
        stop_event.clear()

# Display execution controls
st.header("Pipeline Control")

# Status indicators
col1, col2, col3 = st.columns(3)

with col1:
    st.metric(
        "Current Round",
        f"{st.session_state.pipeline_status['current_round']}/{st.session_state.pipeline_status['total_rounds']}"
    )

with col2:
    st.metric(
        "Progress",
        f"{st.session_state.pipeline_status['progress']:.1f}%"
    )

with col3:
    status = "Running" if st.session_state.pipeline_status["running"] else "Ready"
    st.metric("Status", status)

# Start/Stop buttons
col4, col5 = st.columns(2)

with col4:
    start_disabled = st.session_state.pipeline_status["running"]
    if st.button("Start Pipeline", type="primary", disabled=start_disabled):
        # Reset log position and container
        st.session_state.pipeline_status["last_log_position"] = 0
        st.session_state.log_container = ""
        # Create a copy of the configuration for the thread
        thread_config = dict(st.session_state.pipeline_config)
        thread_status = st.session_state.pipeline_status
        
        # Reset the stop event
        st.session_state.stop_event.clear()
        
        # Create and start the pipeline thread
        st.session_state.pipeline_thread = threading.Thread(
            target=run_pipeline,
            args=(thread_config, thread_status, st.session_state.stop_event),
            daemon=True  # Make thread daemon so it continues when switching pages
        )
        st.session_state.pipeline_thread.start()
        st.rerun()

with col5:
    stop_disabled = not st.session_state.pipeline_status["running"]
    if st.button("Stop Pipeline", type="secondary", disabled=stop_disabled):
        # Set the stop event
        st.session_state.stop_event.set()
        st.session_state.pipeline_status["running"] = False
        st.warning("Stopping pipeline execution... Please wait.")
        st.rerun()

# Progress bar
if st.session_state.pipeline_status["running"]:
    st.progress(st.session_state.pipeline_status["progress"] / 100)

# Log viewer
st.header("Execution Log")

# Read new log content
new_log_content = read_log_file(Path(st.session_state.pipeline_config["out_dir"]) / "logs" / "quick_pipeline.log")
if new_log_content:
    # Append new content to existing log
    st.session_state.log_container += new_log_content
    # Update progress based on log content
    update_progress_from_log(new_log_content)

# Display log with auto-scroll and syntax highlighting
create_auto_scrolling_text_area(st.session_state.log_container)

# Auto-refresh logs while pipeline is running
if st.session_state.pipeline_status["running"]:
    time.sleep(0.5)  # Reduced delay for more responsive updates
    st.rerun()

# Check thread status and update if needed
if st.session_state.pipeline_thread and st.session_state.pipeline_thread.is_alive():
    time.sleep(0.1)  # Small delay to prevent too frequent updates
    st.rerun()
elif st.session_state.pipeline_thread and not st.session_state.pipeline_thread.is_alive():
    # Thread completed, clean up
    st.session_state.pipeline_thread = None
    if "error" in st.session_state.pipeline_status:
        st.error(f"Pipeline failed: {st.session_state.pipeline_status['error']}")
    else:
        st.success("Pipeline completed successfully!")
    st.rerun()

# Display configuration summary
with st.expander("Current Configuration"):
    st.json(st.session_state.pipeline_config) 