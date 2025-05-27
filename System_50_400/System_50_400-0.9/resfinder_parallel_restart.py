from re import M
#import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
#from sympy import Point, Line
from scipy import interpolate
import sys
from gr_wrapper import ckerr_eql2j, ckerr_j2eql, ckerr_minverse, ckerr_minv2omega
import os             
import math
import globalpars
import time


# Defining variables for random number generator
M = globalpars.GLOBALPAR_M
astar = globalpars.GLOBALPAR_astar
mu_outer = globalpars.GLOBALPAR_mu_outer
mu_inner = globalpars.GLOBALPAR_mu_inner
file_amount = globalpars.GLOBALPAR_N_sys
int_min = globalpars.GLOBALPAR_N_min_Fourier
int_max = globalpars.GLOBALPAR_N_max_Fourier
system_label = sys.argv[1]

output_file = f"potential_resonances_{system_label}.txt"

# Define lexicographic check function for modes with -1 multiples
def ispos(x):
    return np.lexsort(np.vstack((x, np.zeros_like(x)))[:,::-1].T)[0] == 1

# Load already computed combinations from output file
completed_combinations = set()
if os.path.exists(output_file):
    with open(output_file, 'r') as f_check:
        for line in f_check:
            parts = line.strip().split()
            if len(parts) >= 14:
                try:
                    n_inner = int(parts[9])
                    k_inner = int(parts[10])
                    n_outer = int(parts[11])
                    k_outer = int(parts[12])
                    combo_key = (n_inner, n_outer, k_inner, k_outer)

                    if not all(isinstance(val, int) for val in combo_key):
                         print(f"Non-integer combo_key found: {combo_key} from line: {line.strip()}")
                    
                    completed_combinations.add(combo_key)
                except ValueError:
                    continue
    print(f"Restarting: {len(completed_combinations)} combinations already processed.")
else:
    print("No existing output file found — starting from scratch.")

# Estimate total number of combinations (including all even/odd filtering)
all_possible_combinations = 0
for n_inner in range(int_min, int_max):
    for n_outer in range(int_min, int_max):
        for k_inner in range(int_min, int_max):
            for k_outer in range(int_min, int_max):
                if abs(k_inner - k_outer) % 2 == 0:
                    all_possible_combinations += 1

remaining = all_possible_combinations - len(completed_combinations)
print(f"Total combinations: {all_possible_combinations}")
print(f"Remaining to process: {remaining}")
print("Starting interpolation...")


# Starting the runs
number = 0
total_checked = 0

# Initialize number from existing output file
if os.path.exists(output_file):
    with open(output_file, 'r') as f:
        for line in reversed(f.readlines()):
            if line.strip() and line.split()[0].isdigit():
                number = int(line.split()[0])
                break
    print(f"Continuing from entry number: {number}")

# Load data
_, time_value_1, J_r_inner_list, J_theta_inner_list, J_phi_inner_list, om_inner_r_list, om_inner_theta_list, om_inner_phi_list, delta_t_list_inner = np.loadtxt("outputs_data/J_evolve_inner_" + str(system_label) + ".txt", unpack=True, skiprows=9)
label, time_value_2, J_r_outer_list, J_theta_outer_list, J_phi_outer_list, om_outer_r_list, om_outer_theta_list, om_outer_phi_list, delta_t_list_outer = np.loadtxt("outputs_data/J_evolve_outer_" + str(system_label) + ".txt", unpack=True, skiprows=9)

# Interpolation of J values
J_inner_r_function = interpolate.CubicSpline(time_value_1, J_r_inner_list)
J_inner_theta_function = interpolate.CubicSpline(time_value_1, J_theta_inner_list)
J_inner_phi_function = interpolate.CubicSpline(time_value_1, J_phi_inner_list)
J_outer_r_function = interpolate.CubicSpline(time_value_2, J_r_outer_list)
J_outer_theta_function = interpolate.CubicSpline(time_value_2, J_theta_outer_list)
J_outer_phi_function = interpolate.CubicSpline(time_value_2, J_phi_outer_list)

# Start time
# Time limits
WALL_TIME_SECONDS = 16 * 60 * 60  # 16 hours
EXIT_BUFFER_SECONDS = 5 * 60     # Exit 5 minutes early
start_time = time.time()

# Time axis

# Define a function that creates nsub subsections within linearly spaced intervals in np.linspace
def subsample(X,nsub):

    L = np.size(X)
    Xsub = np.zeros((L-1)*nsub+1)
    for j in range(L-1):
        Xsub[j*nsub:(j+1)*nsub+1] = np.linspace(X[j],X[j+1],nsub+1)
    return Xsub
# t = np.linspace(time_value_1[0], time_value_1[-1], 1000) # Linearly spaced time grid
n_sub = 10 # Number of sub-divisions within an interval
t = subsample(time_value_1,n_sub) # Subsampled time grid
total_steps = len(t)


# Omega output file
omega_file = f"precomputed_omegas_{system_label}.txt"

# Lists to store omega components
om_inner_r_list = []
om_inner_theta_list = []
om_inner_phi_list = []
om_outer_r_list = []
om_outer_theta_list = []
om_outer_phi_list = []

# Load or compute omega values
if os.path.exists(omega_file):
    with open(omega_file, 'r') as f_in:
        lines = [line for line in f_in if not line.startswith("#")]
    already_done = len(lines)
    print(f"Resuming omega precomputation from line {already_done}/{total_steps}...")

    if already_done > 0:
        data = np.loadtxt(omega_file, comments="#")
        om_inner_r_list = list(data[:, 1])
        om_inner_theta_list = list(data[:, 2])
        om_inner_phi_list = list(data[:, 3])
        om_outer_r_list = list(data[:, 4])
        om_outer_theta_list = list(data[:, 5])
        om_outer_phi_list = list(data[:, 6])
else:
    already_done = 0
    print("Starting fresh omega precomputation...")
    with open(omega_file, 'w') as f_out:
        f_out.write("# i t om_in_r om_in_theta om_in_phi om_out_r om_out_theta om_out_phi\n")

# Compute and append remaining omega values
with open(omega_file, 'a') as f_out:
    for i in range(already_done, total_steps):
        J_inner = [J_inner_r_function(t[i]), J_inner_theta_function(t[i]), J_inner_phi_function(t[i])]
        J_outer = [J_outer_r_function(t[i]), J_outer_theta_function(t[i]), J_outer_phi_function(t[i])]
        _, M_in = ckerr_minverse(J_inner, M, astar)
        _, M_out = ckerr_minverse(J_outer, M, astar)
        omega_in = ckerr_minv2omega(M_in)
        omega_out = ckerr_minv2omega(M_out)

        # Append to in-memory lists
        om_inner_r_list.append(omega_in[0])
        om_inner_theta_list.append(omega_in[1])
        om_inner_phi_list.append(omega_in[2])
        om_outer_r_list.append(omega_out[0])
        om_outer_theta_list.append(omega_out[1])
        om_outer_phi_list.append(omega_out[2])

        # Write to file
        row = [int(i), t[i], *omega_in, *omega_out]
        row_str = " ".join(f"{val:.12e}" if isinstance(val, float) else str(val) for val in row)
        f_out.write(f"{row_str}\n")
        f_out.flush()

        # Optional live output
        print(f"[{i+1}/{total_steps}] {row_str}")

if len(om_inner_r_list) == total_steps:
    print("All omega values precomputed and loaded.")
else:
    raise RuntimeError("Mismatch in omega list lengths — something went wrong.")

# Rebuild Omega_* for compatibility with existing code
Omega_inner = list(zip(om_inner_r_list, om_inner_theta_list, om_inner_phi_list))
Omega_outer = list(zip(om_outer_r_list, om_outer_theta_list, om_outer_phi_list))


# Compute A coefficients
A_inner_r_list = []
A_inner_theta_list = []
A_outer_r_list = []
A_outer_theta_list = []

for i in range(total_steps):
    denom = om_inner_phi_list[i] - om_outer_phi_list[i]
    A_inner_r_list.append(om_inner_r_list[i] / denom)
    A_inner_theta_list.append(om_inner_theta_list[i] / denom)
    A_outer_r_list.append(om_outer_r_list[i] / denom)
    A_outer_theta_list.append(om_outer_theta_list[i] / denom)




# # Calculate all Omega and A values
# for i in range(len(t)):
#     J_inner_array.append([J_inner_r_function(t[i]), J_inner_theta_function(t[i]), J_inner_phi_function(t[i])])
#     J_outer_array.append([J_outer_r_function(t[i]), J_outer_theta_function(t[i]), J_outer_phi_function(t[i])])
#     Omega_inner.append(ckerr_minv2omega(ckerr_minverse(J_inner_array[i], M, astar)[1]))
#     Omega_outer.append(ckerr_minv2omega(ckerr_minverse(J_outer_array[i], M, astar)[1]))

#     om_inner_r_function = Omega_inner[i][0]
#     om_inner_theta_function = Omega_inner[i][1]
#     om_inner_phi_function = Omega_inner[i][2]

#     om_outer_r_function = Omega_outer[i][0]
#     om_outer_theta_function = Omega_outer[i][1]
#     om_outer_phi_function = Omega_outer[i][2]

#     A_inner_r_list.append(om_inner_r_function/(om_inner_phi_function-om_outer_phi_function))
#     A_inner_theta_list.append(om_inner_theta_function/(om_inner_phi_function-om_outer_phi_function))
#     A_outer_r_list.append(om_outer_r_function/(om_inner_phi_function-om_outer_phi_function))
#     A_outer_theta_list.append(om_outer_theta_function/(om_inner_phi_function-om_outer_phi_function))
#     print(f"Loop number: {i}/{len(t)}")

# print("Loop time:", time.time() - start_time)


# Begin resonance search
label.tolist()
[int(num) for num in label]
print("Before nested loops")
with open(output_file, 'a') as f:
    for n_inner in range(int_min, int_max): 
        for n_outer in range(int_min, int_max):
            for k_inner in range(int_min, int_max):
                for k_outer in range(int_min, int_max):

                    # combo_key = (n_inner, n_outer, k_inner, k_outer)
                    # if combo_key in completed_combinations:
                    #     continue

                    combo_key = (n_inner, n_outer, k_inner, k_outer)
                    if combo_key in completed_combinations:
                        print("moving on")
                        continue
                    total_checked += 1
                    print(f"Processing {total_checked} / {remaining} new combinations...", end='\r')

                    if abs(k_inner - k_outer) % 2 == 0:
                        m = []
                        for j in range(len(A_inner_r_list)):
                            m.append(-n_inner*A_inner_r_list[j] - k_inner*A_inner_theta_list[j] + n_outer*A_outer_r_list[j] + k_outer*A_outer_theta_list[j])

                        for i in range(len(m)-1):
                            if math.floor(m[i+1]) != math.floor(m[i]):
                                if abs(m[i]) <= int_max:
                                    m_i = int(np.round((m[i] + m[i+1])/2))
                                    gd = np.gcd.reduce([n_inner, n_outer, k_inner, k_outer, m_i])
                                    if gd < 3 and ispos([m_i, n_inner, n_outer, k_inner, k_outer]):
                                        if(n_inner == 0 and k_inner == 0 and n_outer == 0 and k_outer == 0 and m_i == 0):
                                                continue
                                        k_checker = ((k_inner - k_outer)/2) % 2
                                        if gd == 1 or k_checker != 0:
                                            if t[i] != time_value_1[0]:

                                                omega_inner_after = Omega_inner[i+1][0]*n_inner + Omega_inner[i+1][1]*k_inner + Omega_inner[i+1][2]*(m_i)
                                                omega_inner_before = Omega_inner[i][0]*n_inner + Omega_inner[i][1]*k_inner + Omega_inner[i][2]*(m_i)
                                                delta_t_inner = t[i+1] - t[i]
                                                omega_inner_dot = (omega_inner_after-omega_inner_before)/delta_t_inner

                                                omega_outer_after = Omega_outer[i+1][0]*n_outer + Omega_outer[i+1][1]*k_outer + Omega_outer[i+1][2]*(m_i)
                                                omega_outer_before = Omega_outer[i][0]*n_outer + Omega_outer[i][1]*k_outer + Omega_outer[i][2]*(m_i)
                                                delta_t_outer = t[i+1] - t[i]
                                                omega_outer_dot = (omega_outer_after-omega_outer_before)/delta_t_outer

                                                gamma = -omega_inner_dot + omega_outer_dot

                                                t_res_crossing = t[i] - (-omega_inner_before + omega_outer_before)/gamma

                                                J_inner_res = [J_inner_r_function(t_res_crossing), J_inner_theta_function(t_res_crossing), J_inner_phi_function(t_res_crossing)]
                                                J_outer_res = [J_outer_r_function(t_res_crossing), J_outer_theta_function(t_res_crossing), J_outer_phi_function(t_res_crossing)]

                                                _, EQL_outer = ckerr_j2eql(J_outer_res, M, astar)
                                                _, EQL_inner = ckerr_j2eql(J_inner_res, M, astar)

                                                _, _, anc_outer = ckerr_eql2j(list(EQL_outer), M, astar)
                                                _, _, anc_inner = ckerr_eql2j(list(EQL_inner), M, astar)

                                                _, M_outer = ckerr_minverse(J_outer_res, M, astar)
                                                omega_outer = ckerr_minv2omega(M_outer)
                                                _, M_inner = ckerr_minverse(J_inner_res, M, astar)
                                                omega_inner = ckerr_minv2omega(M_inner)

                                                Delta_omega = n_outer * omega_outer[0] + k_outer * omega_outer[1] + m_i * omega_outer[2] - (n_inner * omega_inner[0] + k_inner * omega_inner[1] + m_i * omega_inner[2])
                                                
                                                number += 1
                                                f.write(f"{number} {system_label} {t[i]} {omega_inner_dot} {omega_outer_dot} {mu_inner} {mu_outer} {gamma} {m_i} {n_inner} {k_inner} {n_outer} {k_outer} {J_inner_res[0]} {J_inner_res[1]} {J_inner_res[2]} {J_outer_res[0]} {J_outer_res[1]} {J_outer_res[2]} {omega_inner[0]} {omega_inner[1]} {omega_inner[2]} {omega_outer[0]} {omega_outer[1]} {omega_outer[2]} {anc_inner[0]} {anc_inner[1]} {anc_inner[2]} {anc_outer[0]} {anc_outer[1]} {anc_outer[2]} {t_res_crossing} {Delta_omega}\n")
                                                f.flush()
                                                os.fsync(f.fileno())
                                                # Check for nearing wall time
                                                elapsed_time = time.time() - start_time
                                                if elapsed_time > (WALL_TIME_SECONDS - EXIT_BUFFER_SECONDS):
                                                    print("\n[INFO] Approaching wall time limit. Exiting cleanly to allow restart.")
                                                    sys.exit(1)
    
                                                # # Optional: write flag file to indicate partial completion
                                                # with open(f"exit_flag_{system_label}.txt", "w") as f:
                                                #     f.write("Clean exit due to time constraint.\n")
                                                #     sys.exit(0)
                                        else:
                                            continue

print("Completed resonance finder")
