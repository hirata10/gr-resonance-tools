import numpy
from subprocess import Popen, PIPE
import sys
import glob, os
from math import ceil
import globalpars
import signal
import time

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

# List to keep track of subprocesses
processes_inner = []
processes_outer = []

# Define the walltime (in seconds) for Chicoma
WALLTIME_LIMIT = 15 * 60 * 60 + 55 * 60  # 15:55 hours, slightly less than walltime for Chicoma at 16 hours (adjust this value)
start_time = time.time()

# Function to handle termination and cleanup
def terminate_processes(signum=None, frame=None):
    print("Walltime exceeded or error occurred, terminating subprocesses...", flush=True)
    # Terminate all subprocesses
    for proc in processes_inner + processes_outer:
        proc.terminate()
    
    # Add a timeout mechanism for process termination
    for proc in processes_inner + processes_outer:
        try:
            proc.wait(timeout=10)  # Wait for 10 seconds
        except subprocess.TimeoutExpired:
            proc.kill()  # Force kill if timeout exceeded
    sys.exit(137)  # Exit the Python script for exceeding walltime

# Set up the signal handler for SIGTERM (sent by SLURM when job is terminated)
signal.signal(signal.SIGTERM, terminate_processes)

# Set up the signal handler for other types of exceptions (e.g., crashes)
def handle_exception(exc_type, exc_value, exc_tb):
    print(f"Uncaught exception: {exc_type}, {exc_value}", flush=True)
    terminate_processes(None, None)
    
# Set the custom exception handler to catch all uncaught exceptions
sys.excepthook = handle_exception


#Defining data for RK4 evolver
t0 = globalpars.GLOBALPAR_t0
n_inner = globalpars.GLOBALPAR_n_time_inner
n_outer = globalpars.GLOBALPAR_n_time_outer
i_start_inner = globalpars.GLOBALPAR_i_start_inner
i_start_outer = globalpars.GLOBALPAR_i_start_outer
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
        
        # Read checkpoint files if they exist, otherwise start from the beginning
        checkpoint_inner_file = f"outputs_data/J_evolve_inner_{int(row[0])}.txt"
        checkpoint_outer_file = f"outputs_data/J_evolve_outer_{int(row[0])}.txt"

        # If the file already exists, then use the last line of data as inputs for the next run and append to that same file
        if os.path.exists(checkpoint_inner_file):
            with open(checkpoint_inner_file, 'r') as f:
                lines_inner = f.readlines()
                # lines_inner = [line.strip() for line in lines_inner if line.strip()]

                if lines_inner:
                    last_line_inner = lines_inner[-1].strip()
                    if (len(last_line_inner.split()) == 9): # Checks that the last line has the right data structure (not plunged)
                        print(f"Last line from the file inner_{int(row[0])}: {last_line_inner}", flush=True)  # Debugging output
                        last_processed_row_inner, last_t_inner, last_J_r_inner, last_J_theta_inner, last_J_phi_inner, last_Omega_r, last_Omega_theta, last_Omega_phi, last_dt = last_line_inner.split()
                        step_inner = n_inner  # Number of total steps
                        which_step = int(last_processed_row_inner) + 1 # Which step to continue from
                        t_inner = float(last_t_inner) # Time to continue from
                        J_r_inner = float(last_J_r_inner)
                        J_theta_inner = float(last_J_theta_inner)
                        J_phi_inner = float(last_J_phi_inner)
                        remaining_steps = step_inner - which_step
                        dt = float(last_dt)
                        label_inner = int(row[0])
                        print(f"The inputs from the restart for inner_{int(row[0])} are: J_r_inner = {J_r_inner}, J_theta_inner = {J_theta_inner}, J_phi_inner = {J_phi_inner}, t_inner = {t_inner}, steps remaining = {remaining_steps}, next step = {which_step}", flush=True)
                        cmd_inner = f"./J_evolve_single_restart {J_r_inner} {J_theta_inner} {J_phi_inner} {t_inner} {step_inner} {globalpars.GLOBALPAR_mu_inner} {label_inner} {'inner'} {which_step} >> outputs_data/J_evolve_inner_{label_inner}.txt"
                        if (which_step < int(n_inner)): # If the number of steps remaining on the restart is less than the total number of steps asked for, then re-run the command
                            processes_inner.append(Popen(cmd_inner, shell = True))
                        else:
                            print(f"End of the requested time steps for inner_{int(row[0])}", flush=True)
                    else:
                        print(f"Last line of data either plunged or is not valid for {checkpoint_inner_file}", flush=True)
                else:
                    print(f"No valid restart data in {checkpoint_inner_file}", flush=True)
                # print(f"Resuming inner body from row {last_processed_row_inner + 1} with last t={last_t_inner}, J_r={last_J_r_inner}, J_theta={last_J_theta_inner}, J_phi={last_J_phi_inner}")
       # If the file does not exist, run using data from action_angle_pairs.txt
        else:
            # Start from the first row if no checkpoint exists
            print(f"There is no initial data, opening the file {checkpoint_inner_file}", flush=True)
            step_inner = n_inner
            which_step = i_start_inner
            t_inner = t0  # Initial time
            J_r_inner = float(row[1]) # {float(row[1])}
            J_theta_inner = float(row[2]) # {float(row[2])}
            J_phi_inner = float(row[3]) # {float(row[3])}
            label_inner = int(row[0])
            cmd_inner = f"./J_evolve_single_norestart {J_r_inner} {J_theta_inner} {J_phi_inner} {t_inner} {step_inner} {globalpars.GLOBALPAR_mu_inner} {label_inner} {'inner'} {which_step} >> outputs_data/J_evolve_inner_{label_inner}.txt"
            processes_inner.append(Popen(cmd_inner, shell = True))
        
        if os.path.exists(checkpoint_outer_file):
            with open(checkpoint_outer_file, 'r') as f:
                lines_outer = f.readlines()
                # lines_outer = [line.strip() for line in lines_outer if line.strip()]

                if lines_outer:
                    last_line_outer = lines_outer[-1].strip()
                    if (len(last_line_outer.split()) == 9): # Checks that the last line has the right data structure (not plunged)
                        print(f"Last line from the file outer_{int(row[0])}: {last_line_outer}", flush=True)
                        last_processed_row_outer, last_t_outer, last_J_r_outer, last_J_theta_outer, last_J_phi_outer, last_Omega_r, last_Omega_theta, last_Omega_phi, last_dt = last_line_outer.split()
                        which_step = int(last_processed_row_outer) + 1
                        step_outer = n_outer
                        t_outer = float(last_t_outer)
                        J_r_outer = float(last_J_r_outer)
                        J_theta_outer = float(last_J_theta_outer)
                        J_phi_outer = float(last_J_phi_outer)
                        remaining_steps = step_outer - which_step
                        dt = float(last_dt)
                        label_outer = int(row[0])
                        print(f"The inputs from the restart for outer_{int(row[0])} are: J_r_outer = {J_r_outer}, J_theta_outer = {J_theta_outer}, J_phi_outer = {J_phi_outer}, t_outer = {t_outer}, steps remaining = {remaining_steps}", flush=True)
                        cmd_outer = f"./J_evolve_single_restart {J_r_outer} {J_theta_outer} {J_phi_outer} {t_outer} {step_outer} {globalpars.GLOBALPAR_mu_outer} {label_outer} {'outer'} {which_step} >> outputs_data/J_evolve_outer_{label_outer}.txt"
                        if (which_step < int(n_outer)):
                            processes_outer.append(Popen(cmd_outer, shell = True))
                        else:
                            print(f"End of the requested time steps for outer_{int(row[0])}", flush=True)
                    else:
                        print(f"Last line of data either plunged or is not valid for {checkpoint_outer_file}", flush=True)
                else:
                    print(f"No valid restart data in {checkpoint_outer_file}", flush=True)
                # print(f"Resuming outer body from row {last_processed_row_outer + 1} with last t={last_t_outer}, J_r={last_J_r_outer}, J_theta={last_J_theta_outer}, J_phi={last_J_phi_outer}")
        else:
            print(f"There is no initial data, opening the file {checkpoint_outer_file}", flush=True)
            step_outer = n_outer
            which_step = i_start_outer
            t_outer = t0  # Initial time
            J_r_outer = float(row[4]) # {float(row[4])} 
            J_theta_outer = float(row[5]) # {float(row[5])} 
            J_phi_outer = float(row[6]) # {float(row[6])}
            label_outer = int(row[0])
            cmd_outer = f"./J_evolve_single_norestart {J_r_outer} {J_theta_outer} {J_phi_outer} {t_outer} {step_outer} {globalpars.GLOBALPAR_mu_outer} {label_outer} {'outer'} {which_step} >> outputs_data/J_evolve_outer_{label_outer}.txt"
            processes_outer.append(Popen(cmd_outer, shell = True))

        # cmd_inner = f"./J_evolve_single {J_r_inner} {J_theta_inner} {J_phi_inner} {t0} {n} {globalpars.GLOBALPAR_mu_inner} {step_inner} {'inner'} >> J_evolve_{"inner"}_{int(row[0])}.txt"
        # cmd_outer = f"./J_evolve_single {J_r_outer} {J_theta_outer} {J_phi_outer} {t0} {n} {globalpars.GLOBALPAR_mu_outer} {step_outer} {'outer'} >> J_evolve_{"outer"}_{int(row[0])}.txt"
        # # fout = open(f"Outer_Body_1/Delta_J_log_{tot_chunk}_{chunk + 1}.txt", "a")
        # processes.append(Popen(cmd, stdout = fout, shell = True))
        # processes_inner.append(Popen(cmd_inner, shell = True))
        # processes_outer.append(Popen(cmd_outer, shell = True))
    
    # Monitor walltime and subprocess completion
    all_processes = processes_inner + processes_outer
    while all_processes:
        # Check if any process is still running
        for process in all_processes[:]:
            if process.poll() is not None:  # If the process has finished, remove it from the list
                all_processes.remove(process)

        # Check the walltime
        elapsed_time = time.time() - start_time
        print(f"Elapsed time: {elapsed_time} seconds", flush=True)  # Debugging
        if elapsed_time > WALLTIME_LIMIT:
            print("Walltime exceeded. Terminating all processes.", flush=True)
            for proc in processes_inner + processes_outer:
                proc.terminate()
            sys.exit(137)  # Exit the script if walltime exceeded

        # Sleep for a bit to avoid checking constantly (e.g., check every 5 seconds)
        time.sleep(100)

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
    
    # Check if walltime is exceeded
    elapsed_time = time.time() - start_time
    if elapsed_time > WALLTIME_LIMIT:
        print("Walltime exceeded. Terminating all processes.", flush=True)
        for proc in processes_inner + processes_outer:
            proc.terminate()
        sys.exit(1)  # Exit the loop and stop further processing

    #fout.close()
sys.exit()

import os
import signal
import subprocess
import sys
import time

# List to keep track of subprocesses
processes_inner = []
processes_outer = []

# Define the walltime (in seconds) (example: 16 hours)
WALLTIME_LIMIT = 16 * 60 * 60  # 16 hours (adjust this value)
start_time = time.time()

# Function to handle termination and cleanup
def terminate_processes(signum=None, frame=None):
    print("Walltime exceeded or error occurred, terminating subprocesses...")
    # Terminate all subprocesses
    for proc in processes_inner + processes_outer:
        proc.terminate()
    sys.exit(0)  # Exit the Python script

# Set up the signal handler for SIGTERM (sent by SLURM when job is terminated)
signal.signal(signal.SIGTERM, terminate_processes)

# Set up the signal handler for other types of exceptions (e.g., crashes)
def handle_exception(exc_type, exc_value, exc_tb):
    print(f"Uncaught exception: {exc_type}, {exc_value}")
    terminate_processes(None, None)
    
# Set the custom exception handler to catch all uncaught exceptions
sys.excepthook = handle_exception

# Executing the compiled code in parallel for both inner and outer body
for chunk in range(0, len(J_data), chunk_size):
    processes_inner = []
    processes_outer = []

    for row in J_data[chunk:chunk + chunk_size]:
        
        # Read checkpoint files if they exist, otherwise start from the beginning
        checkpoint_inner_file = f"outputs_data/J_evolve_inner_{int(row[0])}.txt"
        checkpoint_outer_file = f"outputs_data/J_evolve_outer_{int(row[0])}.txt"

        # Process inner body
        if os.path.exists(checkpoint_inner_file):
            with open(checkpoint_inner_file, 'r') as f:
                lines_inner = f.readlines()

                if lines_inner:
                    last_line_inner = lines_inner[-1].strip()
                    print(f"Last line from the file inner_{int(row[0])}: {last_line_inner}")
                    last_processed_row_inner, last_t_inner, last_J_r_inner, last_J_theta_inner, last_J_phi_inner, last_Omega_r, last_Omega_theta, last_Omega_phi, last_dt = last_line_inner.split()
                    step_inner = n  # Number of total steps
                    which_step = int(last_processed_row_inner) + 1  # Which step to continue from
                    t_inner = float(last_t_inner)  # Time to continue from
                    J_r_inner = float(last_J_r_inner)
                    J_theta_inner = float(last_J_theta_inner)
                    J_phi_inner = float(last_J_phi_inner)
                    remaining_steps = step_inner - which_step
                    dt = float(last_dt)
                    label_inner = int(row[0])
                    print(f"The inputs from the restart for inner_{int(row[0])} are: J_r_inner = {J_r_inner}, J_theta_inner = {J_theta_inner}, J_phi_inner = {J_phi_inner}, t_inner = {t_inner}, steps remaining = {remaining_steps}, next step = {which_step}")
                    cmd_inner = f"./J_evolve_single_restart {J_r_inner} {J_theta_inner} {J_phi_inner} {t_inner} {step_inner} {globalpars.GLOBALPAR_mu_inner} {label_inner} {'inner'} {which_step} >> outputs_data/J_evolve_inner_{label_inner}.txt"
                    if (which_step != int(n) - 1):  # If steps remaining is less than requested steps, run the command
                        processes_inner.append(subprocess.Popen(cmd_inner, shell=True))
                    else:
                        print(f"End of the requested time steps for inner_{int(row[0])}")
                else:
                    print(f"No valid restart data in {checkpoint_inner_file}")
        else:
            # Start from scratch if no restart data exists
            print(f"There is no initial data, opening the file {checkpoint_inner_file}")
            step_inner = n
            which_step = i_start
            t_inner = t0  # Initial time
            J_r_inner = float(row[1])
            J_theta_inner = float(row[2])
            J_phi_inner = float(row[3])
            label_inner = int(row[0])
            cmd_inner = f"./J_evolve_single_norestart {J_r_inner} {J_theta_inner} {J_phi_inner} {t_inner} {step_inner} {globalpars.GLOBALPAR_mu_inner} {label_inner} {'inner'} {which_step} >> outputs_data/J_evolve_inner_{label_inner}.txt"
            processes_inner.append(subprocess.Popen(cmd_inner, shell=True))

        # Process outer body
        if os.path.exists(checkpoint_outer_file):
            with open(checkpoint_outer_file, 'r') as f:
                lines_outer = f.readlines()

                if lines_outer:
                    last_line_outer = lines_outer[-1].strip()
                    print(f"Last line from the file outer_{int(row[0])}: {last_line_outer}")
                    last_processed_row_outer, last_t_outer, last_J_r_outer, last_J_theta_outer, last_J_phi_outer, last_Omega_r, last_Omega_theta, last_Omega_phi, last_dt = last_line_outer.split()
                    which_step = int(last_processed_row_outer) + 1
                    step_outer = n
                    t_outer = float(last_t_outer)
                    J_r_outer = float(last_J_r_outer)
                    J_theta_outer = float(last_J_theta_outer)
                    J_phi_outer = float(last_J_phi_outer)
                    remaining_steps = step_outer - which_step
                    dt = float(last_dt)
                    label_outer = int(row[0])
                    print(f"The inputs from the restart for outer_{int(row[0])}: J_r_outer = {J_r_outer}, J_theta_outer = {J_theta_outer}, J_phi_outer = {J_phi_outer}, t_outer = {t_outer}, steps remaining = {remaining_steps}")
                    cmd_outer = f"./J_evolve_single_restart {J_r_outer} {J_theta_outer} {J_phi_outer} {t_outer} {step_outer} {globalpars.GLOBALPAR_mu_outer} {label_outer} {'outer'} {which_step} >> outputs_data/J_evolve_outer_{label_outer}.txt"
                    if (which_step != int(n) - 1):
                        processes_outer.append(subprocess.Popen(cmd_outer, shell=True))
                    else:
                        print(f"End of the requested time steps for outer_{int(row[0])}")
        else:
            print(f"There is no initial data, opening the file {checkpoint_outer_file}")
            step_outer = n
            which_step = i_start
            t_outer = t0  # Initial time
            J_r_outer = float(row[4])
            J_theta_outer = float(row[5])
            J_phi_outer = float(row[6])
            label_outer = int(row[0])
            cmd_outer = f"./J_evolve_single_norestart {J_r_outer} {J_theta_outer} {J_phi_outer} {t_outer} {step_outer} {globalpars.GLOBALPAR_mu_outer} {label_outer} {'outer'} {which_step} >> outputs_data/J_evolve_outer_{label_outer}.txt"
            processes_outer.append(subprocess.Popen(cmd_outer, shell=True))

    # Wait for all outer processes to finish
    for process in processes_outer:
        process.wait()

    # Check if walltime is exceeded
    elapsed_time = time.time() - start_time
    if elapsed_time > WALLTIME_LIMIT:
        print("Walltime exceeded. Terminating all processes.")
        for proc in processes_inner + processes_outer:
            proc.terminate()
        break  # Exit the loop and stop further processing


sys.exit()

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