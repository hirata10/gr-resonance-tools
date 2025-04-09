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
processes = []
system_number = sys.argv[1] # Which system to compute
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
    sys.exit(137)  # Exit the Python script for exceeding walltime

# Set up the signal handler for SIGTERM (sent by SLURM when job is terminated)
signal.signal(signal.SIGTERM, terminate_processes)

Delta_J_data = numpy.loadtxt(f"Output_Delta_J_{system_number}/final_Delta_J_{system_number}.txt", delimiter = " ")

# Check if data has only one row
if Delta_J_data.ndim == 1:  # If it's a 1D array
    Delta_J_data = Delta_J_data.reshape(1, -1)  # Reshape to 2D with one row

print(len(Delta_J_data))

iCount = 0
chunk_size = globalpars.GLOBALPAR_chunk_size # number of jobs to do in parallel
tot_chunk = ceil(len(Delta_J_data)/chunk_size)

output_dir = f"Output_Delta_Phi_{system_number}"

# Check if the directory exists, and only create it if it does not
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

#For loop that chunks the entire data file into chunk size and run the cmd for the entire chunk size before moving to the next chunk
for chunk in range(0,len(Delta_J_data),chunk_size):
    #print(chunk, chunk + chunk_size)
    
    for row in Delta_J_data[chunk : chunk + chunk_size]:
        cmd = f"./Delta_Phi_single {float(spin)} {float(row[15])} {float(row[16])} {float(row[17])} {float(mass)} {float(row[21])} {float(row[22])} {float(row[23])} {int(row[2])} {int(row[3])} {int(row[6])} {int(row[0])} {int(row[1])}"
        fout = open(f"Output_Delta_Phi_{system_number}/Delta_Phi_{system_number}_log_{tot_chunk}_{chunk + 1}.txt", "a")
        processes.append(Popen(cmd, stdout = fout, shell = True))

        #{J_inner[0]} {J_inner[1]} {J_inner[2]} {J_outer[0]} {J_outer[1]} {J_outer[2]} {omega_inner[0]} {omega_inner[1]} {omega_inner[2]} {omega_outer[0]} {omega_outer[1]} {omega_outer[2]}
    
    for process in processes:
        process.wait()

    fout.close()

# This will concatenate the files from the for loops
with open(f"Output_Delta_Phi_{system_number}/tot_Delta_Phi_{system_number}.txt", 'w') as output_file:
        for filename in os.listdir("Output_Delta_Phi_" + str(system_number) + "/"):
            file_path = os.path.join("Output_Delta_Phi_" + str(system_number) + "/", filename)
            if os.path.isfile(file_path) and filename.endswith('.txt'):
                with open(file_path, 'r') as input_file:
                    output_file.write(input_file.read())
                    output_file.write('\n')  # Add a newline between concatenated files

# This will reorganize the concatenated arrays by the first column
data = numpy.loadtxt("Output_Delta_Phi_" + str(system_number) + "/tot_Delta_Phi_" + str(system_number) + ".txt")

# Check if data is 1D and reshape to 2D (1 row, N columns) if necessary
if data.ndim == 1:  # If it's a 1D array (single row with multiple columns)
    data = data.reshape(1, -1)  # Reshape it to 2D with one row

# Verify the shape of data after reshaping
print("Shape of data after concatenation and reshaping:", data.shape)

sorted_data = data[data[:, 0].argsort()]
numpy.savetxt("Output_Delta_Phi_" + str(system_number) + "/tot_Delta_Phi_" + str(system_number) + ".txt", sorted_data, fmt='%g', delimiter=' ')