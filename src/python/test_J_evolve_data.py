import numpy as np
#import matplotlib
import sys
from gr_wrapper import ckerr_eql2j, ckerr_j2eql, ckerr_minverse, ckerr_minv2omega, j_dot_selfforce
import sys
import globalpars

file = sys.argv[1]

M = globalpars.GLOBALPAR_M
astar = globalpars.GLOBALPAR_astar
nl = globalpars.GLOBALPAR_nl_self
nmax = globalpars.GLOBALPAR_nmax
kmax = globalpars.GLOBALPAR_kmax
mmax = globalpars.GLOBALPAR_mmax

label, time_value_1, J_r_inner_list, J_theta_inner_list, J_phi_inner_list, om_inner_r_list, om_inner_theta_list, om_inner_phi_list, delta_t_list_inner = np.loadtxt(f"outputs_data/{file}.txt", unpack=True, skiprows=9)

length = len(time_value_1)
J_inner = []
EQL_inner = []
anc_inner = []
eccentricity = []
plot_data = []
time = []
inclination = []

for i in range(1):
    t = time_value_1[i]
    J_r_inner = J_r_inner_list[i]
    J_theta_inner = J_theta_inner_list[i]
    J_phi_inner = J_phi_inner_list[i]

    time.append(t)
    J_inner.append([J_r_inner, J_theta_inner, J_phi_inner])

    _, EQL = ckerr_j2eql(J_inner[i], M, astar)
    EQL_inner.append(EQL)
    _, _, anc = ckerr_eql2j(list(EQL_inner[i]), M, astar)

    anc_inner.append(anc)

    eccen = (anc_inner[i][2] - anc_inner[i][1]) / (anc_inner[i][2] + anc_inner[i][1])

    eccentricity.append(eccen)

    # print("Starting J_dotsf")
    # J_dotsf = j_dot_selfforce(nl, nmax, kmax, mmax, anc_inner[i][2], anc_inner[i][1], 0., anc_inner[i][0], M, astar)
    print(time[i], anc_inner[i][0], anc_inner[i][1], anc_inner[i][2], eccentricity[i], flush=True)
    # print(time[i], anc_inner[i][0], anc_inner[i][1], anc_inner[i][2], eccentricity[i], J_dotsf[0], J_dotsf[1], J_dotsf[2], flush=True)


