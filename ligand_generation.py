# ligand_generation.py
import subprocess
import threading
import time
import os
import signal
import resource
import io
import sys

def run_ligand_generation(checkpoint, pdbfile, outfile, resi_list, n_samples, sanitize, log_callback=None):
    # Get the current directory
    current_dir = os.getcwd()
    
    # Normalize paths to be absolute
    checkpoint_path = os.path.abspath(checkpoint)
    pdbfile_path = os.path.abspath(pdbfile)
    outfile_path = os.path.abspath(outfile)
    
    # Change to the directory containing DiffSBDD
    diffsbdd_dir = os.path.join(current_dir, "DiffSBDD")
    
    command = (
        f"cd {diffsbdd_dir} && python generate_ligands.py {checkpoint_path} "
        f"--pdbfile {pdbfile_path} --outfile {outfile_path} --resi_list {' '.join(resi_list)} "
        f"--n_samples {n_samples} {'--sanitize' if sanitize else ''}"
    )
    
    if log_callback:
        log_callback(f"Executing: {command}\nGenerating {n_samples} ligands...\n")

    def run_process():
        start_time = time.time()
        process = None
        stdout_buffer = io.StringIO()
        stderr_buffer = io.StringIO()
        
        try:
            # Try to increase the soft limit for file descriptors for this process
            try:
                soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
                resource.setrlimit(resource.RLIMIT_NOFILE, (min(hard, 4096), hard))
                if log_callback:
                    log_callback(f"Set file descriptor limit to {min(hard, 4096)}")
            except Exception as e:
                if log_callback:
                    log_callback(f"Warning: Could not set file descriptor limit: {e}")
            
            # Use a with statement to ensure proper cleanup
            with subprocess.Popen(
                command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                text=True, close_fds=True, bufsize=1
            ) as process:
                # Create a safe logging function that won't throw exceptions
                def safe_log(message):
                    try:
                        if log_callback and message.strip():
                            log_callback(message.strip())
                    except Exception as e:
                        # If logging fails, write to a local buffer instead
                        print(f"Logging error: {e}", file=sys.stderr)
                        print(f"Message: {message}", file=sys.stderr)
                
                # Read output in chunks to avoid buffer issues
                for line in iter(process.stdout.readline, ''):
                    stdout_buffer.write(line)
                    safe_log(line)
                
                # Read error output
                for err in iter(process.stderr.readline, ''):
                    stderr_buffer.write(err)
                    safe_log(err)
                
                # Wait for process to complete with timeout
                try:
                    process.wait(timeout=3600)  # 1 hour timeout
                except subprocess.TimeoutExpired:
                    safe_log("Process timed out after 1 hour, terminating...")
                    try:
                        process.terminate()
                        process.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        safe_log("Process did not terminate gracefully, killing...")
                        process.kill()
            
            elapsed_time = time.time() - start_time
            if log_callback:
                try:
                    log_callback(f"Finished generating ligands! (Time taken: {elapsed_time:.2f} seconds)\n")
                except Exception:
                    print(f"Finished generating ligands! (Time taken: {elapsed_time:.2f} seconds)", file=sys.stderr)
                    
        except Exception as e:
            if log_callback:
                try:
                    log_callback(f"Error in ligand generation process: {e}")
                except Exception:
                    print(f"Error in ligand generation process: {e}", file=sys.stderr)
            
            # Ensure process is terminated if an exception occurs
            if process and process.poll() is None:
                try:
                    process.terminate()
                    process.wait(timeout=5)
                except:
                    try:
                        process.kill()
                    except:
                        pass
    
    thread = threading.Thread(target=run_process)
    thread.daemon = True  # Make thread a daemon so it doesn't block program exit
    thread.start()
    return thread
