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

res_data = numpy.loadtxt("potential_resonances.txt", delimiter = " ")
print(len(res_data))

#Make sure the Delta_J_log.txt file is removed before each NEW run of this Python script
#Otherwise, it will append the new stuff to the old file

# for fname in glob.glob("Delta_J_log.txt"):
#     fname_bool = input(f"delete file [{fname}] (y/n): ")
#     if fname_bool == "y":
#         os.remove(fname)
#     elif fname_bool == "n":
#         print("appending onto existing file")
#     else:
#         print("invalid inputs")

Plist = []
iCount = 0
chunk_size = globalpars.GLOBALPAR_chunk_size # number of jobs to do in parallel
tot_chunk = ceil(len(res_data)/chunk_size)

#For loop that chunks the entire data file into chunk size and run the cmd for the entire chunk size before moving to the next chunk
for chunk in range(0,len(res_data),chunk_size):
    #print(chunk, chunk + chunk_size)
    processes = []
    
    for row in res_data[chunk : chunk + chunk_size]:
        #print(row)
        theta_res_F = random.uniform(0,2* math.pi)
        cmd = f"./Delta_J_single {int(row[9])} {int(row[11])} {int(row[10])} {int(row[12])} {int(row[8])} {int(row[8])} {float(row[27])} {float(row[26])} {float(row[25])} {float(row[30])} {float(row[29])} {float(row[28])} {float(row[7])} {theta_res_F} {float(row[6])} {int(row[1])} {int(row[0])} {float(row[19])} {float(row[20])} {float(row[21])} {float(row[22])} {float(row[23])} {float(row[24])} {float(row[13])} {float(row[14])} {float(row[15])} {float(row[16])} {float(row[17])} {float(row[18])}"
        fout = open(f"Output_Delta_J/Delta_J_log_{tot_chunk}_{chunk + 1}.txt", "a")
        processes.append(Popen(cmd, stdout = fout, shell = True))

        #{J_inner[0]} {J_inner[1]} {J_inner[2]} {J_outer[0]} {J_outer[1]} {J_outer[2]} {omega_inner[0]} {omega_inner[1]} {omega_inner[2]} {omega_outer[0]} {omega_outer[1]} {omega_outer[2]}
    
    for process in processes:
        process.wait()

    fout.close()

# This will concatenate the files from the for loops
with open("Output_Delta_J/tot_Delta_J.txt", 'w') as output_file:
        for filename in os.listdir("Output_Delta_J/"):
            file_path = os.path.join("Output_Delta_J/", filename)
            if os.path.isfile(file_path) and filename.endswith('.txt'):
                with open(file_path, 'r') as input_file:
                    output_file.write(input_file.read())
                    output_file.write('\n')  # Add a newline between concatenated files

# This will reorganize the concatenated arrays by the first column
data = numpy.loadtxt("Output_Delta_J/tot_Delta_J.txt")
sorted_data = data[data[:, 0].argsort()]
numpy.savetxt("Output_Delta_J/tot_Delta_J.txt", sorted_data, fmt='%g', delimiter=' ')

# for row in res_data[0:50]:
#     #print(row, "First row")
#     cmd = f"./Delta_J_single {int(row[9])} {int(row[11])} {int(row[10])} {int(row[12])} {int(row[8])} {int(row[8])} {float(row[21])} {float(row[20])} {float(row[19])} {float(row[24])} {float(row[23])} {float(row[22])} {float(row[7])} {float(row[6])} {int(row[1])} {int(row[0])}"
#     print(cmd)
#     fout = open("Delta_J_log.txt", "a")
#     Plist += [Popen(cmd, stdout = fout, shell = True)]
#     # iCount += 1
#     # if iCount >= 20:
#     #     Plist[-20].wait()
#     #Popen(cmd, stdout = fout, shell = True)

# sys.exit()

# cmd = []

# for row in res_data[0:3]:
#     #print(row, "First row")
#     cmd.append(f"./Delta_J_single_set {int(row[9])} {int(row[11])} {int(row[10])} {int(row[12])} {int(row[8])} {int(row[8])} {float(row[21])} {float(row[20])} {float(row[19])} {float(row[24])} {float(row[23])} {float(row[22])} {float(row[7])} {float(row[6])} {int(row[1])}")
#     print(cmd)

# P = Popen(cmd, shell = True)
# sys.exit()
# poll = P.poll()

# if poll is None:
#     # p.subprocess is alive
#     print("Running")

# sys.exit()
# for row in res_data[0:2]:
#     #print(row, "First row")
#     cmd = f"./Delta_J_single_set {int(row[9])} {int(row[11])} {int(row[10])} {int(row[12])} {int(row[8])} {int(row[8])} {float(row[21])} {float(row[20])} {float(row[19])} {float(row[24])} {float(row[23])} {float(row[22])} {float(row[7])} {float(row[6])} {int(row[1])}"
#     print(cmd)
#     Popen(cmd, shell = True)