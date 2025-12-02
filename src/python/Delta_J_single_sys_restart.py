import numpy
from subprocess import Popen
import sys
import os
import math
import random
import globalpars
import signal
import time

"""
Input for Delta_J runs. "row" corresponds to the row in the loaded txt file.
So row[i] corresponds the i-th entry in that given row

row[9] = n_res_inner
row[11] = n_res_outer
row[10] = k_res_inner
row[12] = k_res_outer
row[8] = m_res_inner
row[8] = m_res_outer
row[27] = ra_inner
row[26] = rp_inner
row[25] = inc_inner
row[30] = ra_outer
row[29] = rp_outer
row[28] = inc_outer
row[7] = angular_accel
row[6] = mass_outer
row[1] = system_label
row[0] = resonance_label
row[19] = Omega_inner_res_r
row[20] = Omega_inner_res_theta
row[21] = Omega_inner_res_phi
row[22] = Omega_outer_res_r
row[23] = Omega_outer_res_theta
row[24] = Omega_outer_res_phi
row[13] = J_inner_res_r
row[14] = J_inner_res_theta
row[15] = J_inner_res_phi
row[16] = J_outer_res_r
row[17] = J_outer_res_theta
row[18] = J_outer_res_phi

"""

# Graceful shutdown handler
def handle_exit(sig, frame):
    print("Graceful shutdown requested.")
    sys.exit(1)

signal.signal(signal.SIGINT, handle_exit)
signal.signal(signal.SIGTERM, handle_exit)

# Walltime limit (15h 55m in seconds)
start_time = time.time()
walltime_limit_sec = 15 * 3600 + 55 * 60  # 57300 seconds
# walltime_limit_sec = 30  # 30 seconds

# Input system label
system_label = sys.argv[1]
res_data = numpy.loadtxt(f"potential_resonances_{system_label}.txt", delimiter=" ")

# Ensure 2D array
if res_data.ndim == 1:
    res_data = res_data.reshape(1, -1)

print(f"Total resonance rows: {len(res_data)}")

chunk_size = globalpars.GLOBALPAR_chunk_size
# chunk_size = 10
total_chunks = math.ceil(len(res_data) / chunk_size)
output_dir = f"Output_Delta_J_{system_label}"
os.makedirs(output_dir, exist_ok=True)

#Directory to executables
exe_dir = "td_res_exe"

for chunk_index in range(total_chunks):

    # Check walltime BEFORE starting chunk
    elapsed = time.time() - start_time
    # if elapsed > walltime_limit_sec:
    #     print(f"Walltime limit reached ({elapsed/3600:.2f} hours). Exiting gracefully.")
    #     sys.exit(1)
    
    # ./Delta_J_single ... calculation time, 30 minutes = 1800 seconds
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
    chunk_data = res_data[chunk_start:chunk_end]

    chunk_file = os.path.join(
        output_dir,
        f"Delta_J_{system_label}_log_{total_chunks}_{chunk_index + 1}.txt"
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
        theta_res_F = random.uniform(0, 2 * math.pi)
        cmd = f"{exe_dir}/Delta_J_single {int(row[9])} {int(row[11])} {int(row[10])} {int(row[12])} {int(row[8])} {int(row[8])} " \
              f"{float(row[27])} {float(row[26])} {float(row[25])} {float(row[30])} {float(row[29])} {float(row[28])} " \
              f"{float(row[7])} {theta_res_F} {float(row[6])} {int(row[1])} {int(row[0])} {float(row[19])} " \
              f"{float(row[20])} {float(row[21])} {float(row[22])} {float(row[23])} {float(row[24])} " \
              f"{float(row[13])} {float(row[14])} {float(row[15])} {float(row[16])} {float(row[17])} {float(row[18])}"
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

        # Sleep for a bit to avoid checking constantly (e.g., check every 600 seconds = 10 min)
        # time.sleep(600)

    fout.close()

# Combine all chunk files into one sorted output
combined_file = os.path.join(output_dir, f"tot_Delta_J_{system_label}.txt")
with open(combined_file, 'w') as output_file:
    for filename in sorted(os.listdir(output_dir)):
        file_path = os.path.join(output_dir, filename)
        if (os.path.isfile(file_path) and 
            filename.endswith('.txt') and
            filename != os.path.basename(combined_file)):
            with open(file_path, 'r') as input_file:
                output_file.write(input_file.read())
                output_file.write('\n')

# Load combined file, sort by first column, save back
data = numpy.loadtxt(combined_file)
if data.ndim == 1:
    data = data.reshape(1, -1)

print(f"Final combined data shape: {data.shape}")
sorted_data = data[data[:, 0].argsort()]
numpy.savetxt(combined_file, sorted_data, fmt='%g', delimiter=' ')

print(f"Delta_J completed", flush=True)
