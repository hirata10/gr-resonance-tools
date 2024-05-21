import numpy
from subprocess import Popen, PIPE
import sys
import glob, os
from math import ceil
import globalpars

commands = ["ls -l", "which python", "echo Aloha"]
#processes = [Popen(cmd, shell = True) for cmd in commands]

"""
Input for running RK4 evolver for J_dot of the inner and outer body. "row" corresponds to the row in the loaded txt file.
So row[i] corresponds the i-th entry in that given row

row[1] = J_r_inner
row[2] = J_theta_inner
row[3] = J_phi_inner
row[4] = J_r_outer
row[5] = J_theta_outer
row[6] = J_phi_outer
t0 = initial start time
n = number of time steps


"""

#Defining data for RK4 evolver
t0 = globalpars.GLOBALPAR_t0
n = globalpars.GLOBALPAR_n_time
J_data = numpy.loadtxt("action_angle_pairs.txt", delimiter = " ")
print(len(J_data))

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

#Defining the number of jobs to run in parallel
chunk_size = globalpars.GLOBALPAR_chunk_size # number of jobs to do in parallel
tot_chunk = ceil(len(J_data)/chunk_size)

#For loop that chunks the entire data file into chunk size and run the cmd for the entire chunk size before moving to the next chunk

#Executing the compiled code in parallel for both inner and outer body
for chunk in range(0,len(J_data),chunk_size):
    #print(chunk, chunk + chunk_size)
    processes_inner = []
    processes_outer = []

    for row in J_data[chunk : chunk + chunk_size]:
        #print(row)
        cmd_inner = f"./J_evolve_single {float(row[1])} {float(row[2])} {float(row[3])} {t0} {n} {globalpars.GLOBALPAR_mu_inner} {int(row[0])} {'inner'}"
        cmd_outer = f"./J_evolve_single {float(row[4])} {float(row[5])} {float(row[6])} {t0} {n} {globalpars.GLOBALPAR_mu_outer} {int(row[0])} {'outer'}"
        #fout = open(f"Outer_Body_1/Delta_J_log_{tot_chunk}_{chunk + 1}.txt", "a")
        #processes.append(Popen(cmd, stdout = fout, shell = True))
        processes_inner.append(Popen(cmd_inner, shell = True))
        processes_outer.append(Popen(cmd_outer, shell = True))
    
    for process in processes_outer:
        process.wait()

    #fout.close()
sys.exit()
#Outer body evolution
for chunk in range(0,len(J_data),chunk_size):
    #print(chunk, chunk + chunk_size)
    processes = []

    for row in J_data[chunk : chunk + chunk_size]:
        #print(row)
        cmd_outer = f"./J_evolve {float(row[4])} {float(row[5])} {float(row[6])} {t0} {n}"
        #fout = open(f"Output_Delta_J/Delta_J_log_{tot_chunk}_{chunk + 1}.txt", "a")
        #processes.append(Popen(cmd, stdout = fout, shell = True))
        processes.append(Popen(cmd_outer, shell = True))
    
    for process in processes:
        process.wait()

    #fout.close()

# This will concatenate the files from the for loops
# with open("Output_Delta_J/tot_Delta_J.txt", 'w') as output_file:
#         for filename in os.listdir("Output_Delta_J/"):
#             file_path = os.path.join("Output_Delta_J/", filename)
#             if os.path.isfile(file_path) and filename.endswith('.txt'):
#                 with open(file_path, 'r') as input_file:
#                     output_file.write(input_file.read())
#                     output_file.write('\n')  # Add a newline between concatenated files

# # This will reorganize the concatenated arrays by the first column
# data = numpy.loadtxt("Output_Delta_J/tot_Delta_J.txt")
# sorted_data = data[data[:, 0].argsort()]
# numpy.savetxt("Output_Delta_J/tot_Delta_J.txt", sorted_data, fmt='%g', delimiter=' ')

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