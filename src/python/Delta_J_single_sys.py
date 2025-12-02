import numpy
from subprocess import Popen, PIPE
import sys
import glob, os
from math import ceil
import math
import random
import globalpars

commands = ["ls -l", "which python", "echo Aloha"]
#processes = [Popen(cmd, shell = True) for cmd in commands]

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
system_label = sys.argv[1]

res_data = numpy.loadtxt("potential_resonances_" + str(system_label) + ".txt", delimiter = " ")

# Check if data has only one row
if res_data.ndim == 1:  # If it's a 1D array
    res_data = res_data.reshape(1, -1)  # Reshape to 2D with one row

print(len(res_data))

#Directory to executables
exe_dir = "td_res_exe"

Plist = []
iCount = 0
chunk_size = globalpars.GLOBALPAR_chunk_size # number of jobs to do in parallel
tot_chunk = ceil(len(res_data)/chunk_size)

os.mkdir("Output_Delta_J_" + str(system_label))

#For loop that chunks the entire data file into chunk size and run the cmd for the entire chunk size before moving to the next chunk
for chunk in range(0,len(res_data),chunk_size):
    #print(chunk, chunk + chunk_size)
    processes = []
    
    for row in res_data[chunk : chunk + chunk_size]:
        #print(row)
        theta_res_F = random.uniform(0,2 * math.pi)
        cmd = f"{exe_dir}/Delta_J_single {int(row[9])} {int(row[11])} {int(row[10])} {int(row[12])} {int(row[8])} {int(row[8])} {float(row[27])} {float(row[26])} {float(row[25])} {float(row[30])} {float(row[29])} {float(row[28])} {float(row[7])} {theta_res_F} {float(row[6])} {int(row[1])} {int(row[0])} {float(row[19])} {float(row[20])} {float(row[21])} {float(row[22])} {float(row[23])} {float(row[24])} {float(row[13])} {float(row[14])} {float(row[15])} {float(row[16])} {float(row[17])} {float(row[18])}"
        fout = open(f"Output_Delta_J_{system_label}/Delta_J_{system_label}_log_{tot_chunk}_{chunk + 1}.txt", "a")
        processes.append(Popen(cmd, stdout = fout, shell = True))

        #{J_inner[0]} {J_inner[1]} {J_inner[2]} {J_outer[0]} {J_outer[1]} {J_outer[2]} {omega_inner[0]} {omega_inner[1]} {omega_inner[2]} {omega_outer[0]} {omega_outer[1]} {omega_outer[2]}
    
    for process in processes:
        process.wait()

    fout.close()

# This will concatenate the files from the for loops
with open("Output_Delta_J_" + str(system_label) + "/tot_Delta_J_" + str(system_label) + ".txt", 'w') as output_file:
        for filename in os.listdir("Output_Delta_J_" + str(system_label) + "/"):
            file_path = os.path.join("Output_Delta_J_" + str(system_label) + "/", filename)
            if os.path.isfile(file_path) and filename.endswith('.txt'):
                with open(file_path, 'r') as input_file:
                    output_file.write(input_file.read())
                    output_file.write('\n')  # Add a newline between concatenated files

# This will reorganize the concatenated arrays by the first column
# This will reorganize the concatenated arrays by the first column
data = numpy.loadtxt("Output_Delta_J_" + str(system_label) + "/tot_Delta_J_" + str(system_label) + ".txt")

# Check if data is 1D and reshape to 2D (1 row, N columns) if necessary
if data.ndim == 1:  # If it's a 1D array (single row with multiple columns)
    data = data.reshape(1, -1)  # Reshape it to 2D with one row

# Verify the shape of data after reshaping
print("Shape of data after concatenation and reshaping:", data.shape)

sorted_data = data[data[:, 0].argsort()]
numpy.savetxt("Output_Delta_J_" + str(system_label) + "/tot_Delta_J_" + str(system_label) + ".txt", sorted_data, fmt='%g', delimiter=' ')
