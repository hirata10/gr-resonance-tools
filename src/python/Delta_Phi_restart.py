import numpy
import subprocess
import signal
from subprocess import Popen, PIPE
import sys
import glob, os
from math import ceil
import math
import random
import globalpars
import time

"""
Input for Delta_Phi runs. "row" corresponds to the row in the loaded txt file.
So row[i] corresponds the i-th entry in that given row

spin = spin of SMBH
row[15] = J_res_inner_r
row[16] = J_res_inner_theta
row[17] = J_res_inner_phi
mass = mass of SMBH
row[21] = Delta_J_tidal_r
row[22] = Delta_J_tidal_theta
row[23] = Delta_J_tidal_phi
row[2] = n_inner
row[3] = k_inner
row[6] = m_inner
row[0] = res_label
row[1] = system_label

"""
system_label = sys.argv[1] # Which system to compute
spin = globalpars.GLOBALPAR_astar
mass = globalpars.GLOBALPAR_M

# Function to handle termination and cleanup
def terminate_processes(signum=None, frame=None):
    print("Walltime exceeded or error occurred, terminating subprocesses...", flush=True)
    # Terminate all subprocesses
    for proc in processes:
        proc.terminate()
    
    # Add a timeout mechanism for process termination
    for proc in processes:
        try:
            proc.wait(timeout=10)  # Wait for 10 seconds
        except subprocess.TimeoutExpired:
            proc.kill()  # Force kill if timeout exceeded
    sys.exit(1)  # Exit the Python script for exceeding walltime

# Set up the signal handler for SIGTERM (sent by SLURM when job is terminated)
signal.signal(signal.SIGTERM, terminate_processes)

# Walltime limit (15h 55m in seconds)
start_time = time.time()
walltime_limit_sec = 15 * 3600 + 55 * 60  # 57300 seconds

Delta_J_data = numpy.loadtxt(f"Output_Delta_J_{system_label}/tot_Delta_J_{system_label}.txt", delimiter = " ")

# Check if data has only one row
if Delta_J_data.ndim == 1:  # If it's a 1D array
    Delta_J_data = Delta_J_data.reshape(1, -1)  # Reshape to 2D with one row

print(f"Length of Delta_J_data array: {len(Delta_J_data)}")

chunk_size = globalpars.GLOBALPAR_chunk_size
# chunk_size = 10
total_chunks = math.ceil(len(Delta_J_data) / chunk_size)
output_dir = f"Output_Delta_Phi_{system_label}"
os.makedirs(output_dir, exist_ok=True)

#Directory to executables
exe_dir = "td_res_exe"

for chunk_index in range(total_chunks):

    # Check walltime BEFORE starting chunk
    elapsed = time.time() - start_time
    # if elapsed > walltime_limit_sec:
    #     print(f"Walltime limit reached ({elapsed/3600:.2f} hours). Exiting gracefully.")
    #     sys.exit(1)
    
    # ./Delta_Phi_single ... calculation time, 30 minutes = 1800 seconds
    estimated_chunk_runtime = 1800

    # Check walltime BEFORE starting chunk
    if elapsed > walltime_limit_sec:
        print(f"Walltime limit reached BEFORE chunk {chunk_index + 1}. Exiting.")
        sys.exit(137)
    # Check if there is enough time to compute the next chunk
    if walltime_limit_sec - elapsed < estimated_chunk_runtime:
        print(f"Not enough time left to safely start chunk {chunk_index + 1}. "
              f"Remaining: {(walltime_limit_sec - elapsed)/3600:.1f} sec "
              f"Required: {estimated_chunk_runtime/3600:.1f} sec")
        sys.exit(137)

    chunk_start = chunk_index * chunk_size
    chunk_end = chunk_start + chunk_size
    chunk_data = Delta_J_data[chunk_start:chunk_end]

    chunk_file = os.path.join(
        output_dir,
        f"Delta_Phi_{system_label}_log_{total_chunks}_{chunk_index + 1}.txt"
    )

    expected_lines = len(chunk_data)  # chunk_size or fewer for last chunk

    # Check if chunk output file exists and is complete
    if os.path.exists(chunk_file):
        try:
            if os.path.getsize(chunk_file) == 0:
                # Empty file: assume incomplete
                print(f"Chunk {chunk_index + 1} is empty, rerunning")
                open(chunk_file, "w").close()  # Optional, already empty
            else:
                data = numpy.loadtxt(chunk_file)
                if data.size == 0:
                    print(f"Chunk {chunk_index + 1} contains no data, rerunning")
                    open(chunk_file, "w").close()
                else:
                    if data.ndim == 1:
                        line_count = 1
                    else:
                        line_count = data.shape[0]

                    if line_count >= expected_lines:
                        print(f"Skipping chunk {chunk_index + 1} ({line_count}/{expected_lines} lines)")
                        continue
                    else:
                        print(f"Chunk {chunk_index + 1} incomplete ({line_count}/{expected_lines}), clearing and rerunning")
                        open(chunk_file, "w").close()
        except Exception as e:
            print(f"Error reading {chunk_file}: {e} — clearing and rerunning")
            open(chunk_file, "w").close()


    print(f"Processing chunk {chunk_index + 1} (rows {chunk_start} to {chunk_end - 1})")

    fout = open(chunk_file, "a")
    processes = []

    for row in chunk_data:
        cmd = f"{exe_dir}/Delta_Phi_single {float(spin)} {float(row[15])} {float(row[16])} {float(row[17])} {float(mass)} {float(row[21])} {float(row[22])} {float(row[23])} {int(row[2])} {int(row[3])} {int(row[6])} {int(row[0])} {int(row[1])}"
        processes.append(Popen(cmd, stdout=fout, shell=True))

    # Monitor walltime and subprocess completion
    while processes:
        # Check if any process is still running
        for process in processes[:]:
            if process.poll() is not None:  # If the process has finished, remove it from the list
                processes.remove(process)

        # Check the walltime
        elapsed = time.time() - start_time
        #print(f"Elapsed time: {elapsed} seconds", flush=True)  # Debugging
        if elapsed > walltime_limit_sec:
            print(f"Walltime exceeded limit ({elapsed/3600:.2f} hours). Terminating all processes.", flush=True)
            for proc in processes:
                proc.terminate()
            sys.exit(137)  # Exit the script if walltime exceeded

    fout.close()

# This will concatenate the files from the for loops
with open(f"Output_Delta_Phi_{system_label}/tot_Delta_Phi_{system_label}.txt", 'w') as output_file:
        for filename in os.listdir("Output_Delta_Phi_" + str(system_label) + "/"):
            file_path = os.path.join("Output_Delta_Phi_" + str(system_label) + "/", filename)
            if os.path.isfile(file_path) and filename.endswith('.txt'):
                with open(file_path, 'r') as input_file:
                    output_file.write(input_file.read())
                    output_file.write('\n')  # Add a newline between concatenated files

# This will reorganize the concatenated arrays by the first column
data = numpy.loadtxt("Output_Delta_Phi_" + str(system_label) + "/tot_Delta_Phi_" + str(system_label) + ".txt")

# Check if data is 1D and reshape to 2D (1 row, N columns) if necessary
if data.ndim == 1:  # If it's a 1D array (single row with multiple columns)
    data = data.reshape(1, -1)  # Reshape it to 2D with one row

# Verify the shape of data after reshaping
print("Shape of data after concatenation and reshaping:", data.shape)

sorted_data = data[data[:, 0].argsort()]
numpy.savetxt("Output_Delta_Phi_" + str(system_label) + "/tot_Delta_Phi_" + str(system_label) + ".txt", sorted_data, fmt='%g', delimiter=' ')

print("Delta_Phi completed")