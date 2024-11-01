import numpy
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
system_label = sys.argv[1]
spin = sys.argv[2]
mass = globalpars.GLOBALPAR_M

res_data = numpy.loadtxt("tot_Delta_J/tot_Delta_J_" + str(system_label) + ".txt", delimiter = " ")
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

os.mkdir("Output_Delta_Phi_" + str(system_label))

#For loop that chunks the entire data file into chunk size and run the cmd for the entire chunk size before moving to the next chunk
# for chunk in range(0,len(res_data),chunk_size):
#     #print(chunk, chunk + chunk_size)
#     processes = []
    
#     for row in res_data[chunk : chunk + chunk_size]:
#         #print(row)
#         cmd = f"./Delta_Phi_single  {float(spin)} {float(row[15])} {float(row[16])} {float(row[17])} {float(mass)} {float(row[21])} {float(row[22])} {float(row[23])} {int(row[2])} {int(row[3])} {int(row[6])} {int(row[0])} {int(row[1])}"
#         fout = open(f"Output_Delta_Phi_{system_label}/Delta_Phi_{system_label}_log_{tot_chunk}_{chunk + 1}.txt", "a")
#         processes.append(Popen(cmd, stdout = fout, shell = True))

#         #{J_inner[0]} {J_inner[1]} {J_inner[2]} {J_outer[0]} {J_outer[1]} {J_outer[2]} {omega_inner[0]} {omega_inner[1]} {omega_inner[2]} {omega_outer[0]} {omega_outer[1]} {omega_outer[2]}
    
#     for process in processes:
#         process.wait()

#     fout.close()

for chunk in range(0, len(res_data), chunk_size):
    processes = []
    
    for row_index, row in enumerate(res_data[chunk : chunk + chunk_size]):
        # Create a unique filename for each process's output
        output_filename = f"Output_Delta_Phi_{system_label}/Delta_Phi_{system_label}_log_{tot_chunk}_{chunk + row_index + 1}.txt"
        
        cmd = f"./Delta_Phi_single {float(spin)} {float(row[15])} {float(row[16])} {float(row[17])} {float(mass)} {float(row[21])} {float(row[22])} {float(row[23])} {int(row[2])} {int(row[3])} {int(row[6])} {int(row[0])} {int(row[1])}"
        
        # Open a unique output file for each process
        with open(output_filename, "w") as fout:
            processes.append(Popen(cmd, stdout=fout, shell=True))

    for process in processes:
        process.wait()


# This will concatenate the files from the for loops
with open("Output_Delta_Phi_" + str(system_label) + "/tot_Delta_Phi_" + str(system_label) + ".txt", 'w') as output_file:
        for filename in os.listdir("Output_Delta_Phi_" + str(system_label) + "/"):
            file_path = os.path.join("Output_Delta_Phi_" + str(system_label) + "/", filename)
            if os.path.isfile(file_path) and filename.endswith('.txt'):
                with open(file_path, 'r') as input_file:
                    output_file.write(input_file.read())
                    output_file.write('\n')  # Add a newline between concatenated files

# This will reorganize the concatenated arrays by the first column
data = numpy.loadtxt("Output_Delta_Phi_" + str(system_label) + "/tot_Delta_Phi_" + str(system_label) + ".txt")
sorted_data = data[data[:, 0].argsort()]
numpy.savetxt("Output_Delta_Phi_" + str(system_label) + "/tot_Delta_Phi_" + str(system_label) + ".txt", sorted_data, fmt='%g', delimiter=' ')
