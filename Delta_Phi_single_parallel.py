import concurrent.futures
import multiprocessing as mp
from pathlib import Path
import subprocess
import queue

import click  # Make sure to `pip install click`
import pandas as pd  # Make sure to `pip install pandas`
from tqdm import tqdm  # Make sure to `pip install tqdm`

import globalpars as gp

mp.set_start_method("fork") # FOR MAC OS
# Define Constants
SYSTEM_LABEL: str = gp.GLOBALPAR_system_label
SPIN: float = gp.GLOBALPAR_astar
MASS: float = gp.GLOBALPAR_M

BASE_DIR = Path("./")
EXECUTABLE = BASE_DIR / "Delta_Phi_single"
CONFIG_FILE = BASE_DIR / "tot_Delta_J_{system_label}.txt"
LOG_FILE = BASE_DIR / "progress.tsv"

log_queue = mp.Queue()
stop_logging_event = mp.Event()

def logger_process(log_queue:mp.Queue=log_queue, stop_event:mp.Event=stop_logging_event):
    """A dedicated process to write log entries to the file in a process-safe manner."""
    while not stop_event.is_set():  # Check if the stop signal has been set
        try:
            log_entry = log_queue.get(timeout=1)  # Wait for an item for 1 second
            if log_entry == "STOP":
                break  # Stop the process safely when receiving the "STOP" signal
            # Write the log entry to the log file
            with open(LOG_FILE, "a") as log_file:
                log_file.write(log_entry)
        except queue.Empty:
            continue  # If the queue is empty, just continue checking for stop signal

    print("Logger process has finished.")


def load_configurations(system_label:str) -> pd.DataFrame:
    """Load numerical configurations from the file using Pandas.

    Args:
        system_label (str): system label name for config file.

    Returns:
        pd.DataFrame: Loaded configurations as dataframe with subset of columns.
    """
    filename = str(CONFIG_FILE).format(system_label=system_label)
    df = pd.read_csv(filename, sep="\\s+", header=None, usecols=[0, 1, 2, 3, 6, 15, 16, 17, 21, 22, 23])
    #  Drop records where entire row is empty
    df = df.dropna(how="all").reset_index().drop(columns=["index"])

    # Rename the columns accordingly
    df.columns = [
        "res_label", "system_label", "n_inner", "k_inner", "m_inner",
        "J_res_inner_r", "J_res_inner_theta", "J_res_inner_phi",
        "Delta_J_tidal_r", "Delta_J_tidal_theta", "Delta_J_tidal_phi"
    ]
    return df
    

def load_completed() -> list:
    """Load the completed configurations from the log file.

    Returns:
        list: List of completed indicies
    """
    try:
        # Read the log file as a DataFrame
        log_df = pd.read_csv(LOG_FILE, sep='\t', header=None, index_col=0)
        completed_indices = log_df.index.tolist()
        return completed_indices
    
    except FileNotFoundError:
        # Log the missing log file warning
        print(f"Warning: {LOG_FILE} not found. Returning empty list.")
        return []

def run_script(index:int, spin:float, J_res_inner_r:float, J_res_inner_theta:float, J_res_inner_phi:float, mass:float,
               Delta_J_tidal_r:float, Delta_J_tidal_theta:float, Delta_J_tidal_phi:float,
               n_inner:int, k_inner:int, m_inner:int, res_label:int, system_label:int, smoke_test:bool=False) -> str:
    """Run the C script with the given arguments, log the success and capture the output.

    Args:
        index (int): Configuration number used for tracking progress.
        spin (float): _description_
        J_res_inner_r (float): _description_
        J_res_inner_theta (float): _description_
        J_res_inner_phi (float): _description_
        mass (float): _description_
        Delta_J_tidal_r (float): _description_
        Delta_J_tidal_theta (float): _description_
        Delta_J_tidal_phi (float): _description_
        n_inner (int): _description_
        k_inner (int): _description_
        m_inner (int): _description_
        res_label (int): _description_
        system_label (int): _description_
        smoke_test (bool, optional): Perform smoke test (not actual executable). Defaults to False.

    Returns:
        str: simple log statement.
    """
    
    # Prepare the command
    # In the case of smoke test, override EXECUTABLE
    if smoke_test:
        # Sleep for one second then print
        cmd = "sleep 1 && echo 1 \t 1 \t 1 \t 1 \t 1 \t 1.0 \t 1.0 \t 1.0 \t 1.0 \t 1.0 \t 1.0 \t 1.0 \t 1.0 \t 1.0 \t 1.0"
    else:
        cmd = f"{EXECUTABLE.resolve()} {spin} {J_res_inner_r} {J_res_inner_theta} {J_res_inner_phi} {mass} {Delta_J_tidal_r} {Delta_J_tidal_theta} {Delta_J_tidal_phi} {n_inner} {k_inner} {m_inner} {res_label} {system_label}"
    
    try:
        # Run the command and capture the output
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Capture the stdout output (which contains the 15 values)
        output = result.stdout.strip()

        # Log the result with the index
        log_output = f"{index}\t{output}\n"
        
        # Append the log to the log file
        log_queue.put(log_output)

        # Indicate successful completion
        return f"Index {index}: Success"
    
    except subprocess.CalledProcessError as e:
        # Handle any errors during the script execution
        error_message = e.stderr.strip() if e.stderr else "Unknown error"
        return f"Index {index}: Error - {error_message}"

@click.command()
@click.option("--system-label", type=str, default=str(SYSTEM_LABEL), show_default=True)
@click.option("--spin", type=float, default=float(SPIN), show_default=True)
@click.option("--mass", type=float, default=float(MASS), show_default=True)
@click.option("--cores", type=int, default=-1, show_default=True, help="Number of cores to use. -1 will use all available cores.")
@click.option("--clean", type=bool, default=False, show_default=True, help="Clean run.")
@click.option("--smoke", type=bool, default=False, show_default=True, help="Run smoke test.")
def main(system_label:str, spin:float, mass:float, cores:int|None=None, clean:bool=False, smoke:bool=False):
    """Main entrypoint with optional CLI inputs.

    Args:
        system_label (str): Optional override of system label used to fetch configuration file.
        spin (float): Optional override of spin variable.
        mass (float): Optional override of mass variable.
        cores (int | None, optional): _description_. Defaults to None.
        clean (bool, optional): Whether or not to perform clean run and not load completed. Defaults to False.
        smoke_test (bool, optional): Perform smoke test (not actual executable). Defaults to False.
    """
    # Load configurations from txt file
    configurations = load_configurations(system_label)

    # Fetch completed configurations/results
    # When running in clean mode, set completed to []
    if not clean:
        completed = load_completed()
    else:
        # TODO: Also clear progress file
        completed = []

    # Filter out configurations that are already completed
    pending_configs = configurations.loc[~configurations.index.isin(completed)]

    print(f"Resuming execution with {len(pending_configs)} pending configurations...")

    # Start logging process
    log_process = mp.Process(target=logger_process, args=(log_queue, stop_logging_event))
    log_process.start()

    # Get/set number of workers
    # Remove one worker to save room for logger process
    if cores == -1:
        num_workers = min(mp.cpu_count(), len(pending_configs)) - 1
    else:
        num_workers = max(1, cores - 1)
    
    print("Number of cores available:", num_workers)
    
    # Starting progress bar
    with tqdm(total=len(pending_configs)) as pbar:
        # Run pending configurations in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers, mp_context=mp.get_context("fork")) as executor:
            futures = {}

            for idx, row in pending_configs.iterrows():
                # Construct the dictionary of arguments from the DataFrame row
                # TODO: Add data type casting to load_configuration as dtype
                args_dict = {
                    "index": idx,
                    "spin": spin,
                    "J_res_inner_r": float(row["J_res_inner_r"]),
                    "J_res_inner_theta": float(row["J_res_inner_theta"]),
                    "J_res_inner_phi": float(row["J_res_inner_phi"]),
                    "mass": mass,
                    "Delta_J_tidal_r": float(row["Delta_J_tidal_r"]),
                    "Delta_J_tidal_theta": float(row["Delta_J_tidal_theta"]),
                    "Delta_J_tidal_phi": float(row["Delta_J_tidal_phi"]),
                    "n_inner": int(row["n_inner"]),
                    "k_inner": int(row["k_inner"]),
                    "m_inner": int(row["m_inner"]),
                    "res_label": int(row["res_label"]),
                    "system_label": str(row["system_label"]),
                    "smoke_test": smoke
                }

                # Submit the task to the executor, unpacking the dictionary into **kwargs
                futures[executor.submit(run_script, **args_dict)] = idx

            for _ in concurrent.futures.as_completed(futures):
                # Increment progress bar
                pbar.update(1)
                # print(_.result())

    # Stop logging process safely
    stop_logging_event.set()  # Set the event to notify the logger process to stop
    log_process.join()  # Wait for the logger thread to finish

if __name__ == "__main__":
    main()
