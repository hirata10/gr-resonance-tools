
from re import M
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from sympy import Point, Line
from scipy import interpolate
import sys
from gr_wrapper import ckerr_eql2j, ckerr_j2eql, ckerr_minverse, ckerr_minv2omega
import os             
import math
import globalpars





"""
Hello!

After you simulate all the orbits and put them in a file,
we use this code to generate a list of all potential resonances
for delta_J to determine the strength of.

Inputs: mu_inner, mu_outer, M, astar, file counter

Outputs (by column): resonance counter, file counter, 
resonance time, omega_inner_dot, omega_outer_dot, 
mu_inner, mu_outer, gamma, m, n_inner, k_inner, n_outer, 
k_outer, J_r_inner, J_theta_inner, J_phi_inner, J_r_outer, 
J_theta_outer, J_phi_outer, inner inclination, inner periapse, 
inner apoapse, outer inclination, outer periapse, outer apoapse
"""

#Defining variables for random number generator
M = globalpars.GLOBALPAR_M
astar = globalpars.GLOBALPAR_astar
mu_outer = globalpars.GLOBALPAR_mu_outer
mu_inner = globalpars.GLOBALPAR_mu_inner
file_amount = globalpars.GLOBALPAR_N_sys
int_min = globalpars.GLOBALPAR_N_min_Fourier
int_max = globalpars.GLOBALPAR_N_max_Fourier

# Open a text file for writing
with open('potential_resonances.txt', 'w') as f:

    #Starting the runs
    number = 0
    for file_number in range(file_amount):

        #Importing all the info from the files. Feel free to change this depending on the files you made
        #Also note how we use label, because there will (very likely) always be less data points in outer than inner.
        _, time_value_1, J_r_inner_list, J_theta_inner_list, J_phi_inner_list, om_inner_r_list, om_inner_theta_list, om_inner_phi_list, delta_t_list_inner = np.loadtxt("outputs_data/J_evolve_inner_" + str(file_number+1) + ".txt", unpack=True, skiprows=2)
        label, time_value_2, J_r_outer_list, J_theta_outer_list, J_phi_outer_list, om_outer_r_list, om_outer_theta_list, om_outer_phi_list, delta_t_list_outer = np.loadtxt("outputs_data/J_evolve_outer_" + str(file_number+1) + ".txt", unpack=True, skiprows=2) 
        
        #Interpolation of J values
        J_inner_r_function = interpolate.interp1d(time_value_1, J_r_inner_list, kind='linear')
        J_inner_theta_function = interpolate.interp1d(time_value_1, J_theta_inner_list, kind='linear')
        J_inner_phi_function = interpolate.interp1d(time_value_1, J_phi_inner_list, kind='linear')
        J_outer_r_function = interpolate.interp1d(time_value_2, J_r_outer_list, kind='linear')
        J_outer_theta_function = interpolate.interp1d(time_value_2, J_theta_outer_list, kind='linear')
        J_outer_phi_function = interpolate.interp1d(time_value_2, J_phi_outer_list, kind='linear')

        #Interpolation of omega values
        om_inner_r_function = interpolate.interp1d(time_value_1, om_inner_r_list, kind='linear')
        om_inner_theta_function = interpolate.interp1d(time_value_1, om_inner_theta_list, kind='linear')
        om_inner_phi_function = interpolate.interp1d(time_value_1, om_inner_phi_list, kind='linear')
        om_outer_r_function = interpolate.interp1d(time_value_2, om_outer_r_list, kind='linear')
        om_outer_theta_function = interpolate.interp1d(time_value_2, om_outer_theta_list, kind='linear')
        om_outer_phi_function = interpolate.interp1d(time_value_2, om_outer_phi_list, kind='linear')

        #Setting up time axis. Note that since the inner body will (usually) run to less than the outer body, so we use its timescale.
        #Also, if any inner files goes to less than 100000, feel free to adjust this accordingly.
        t = np.linspace(1, time_value_1[-1], 100000)

        #Composing A values (which are)
        A_inner_r_list = []
        A_inner_theta_list = []
        A_outer_r_list = []
        A_outer_theta_list = []

        res_timing = []

        i = 0
        for i in range(len(t)):

            A_inner_r_list.append(om_inner_r_function(t[i])/(om_inner_phi_function(t[i])-om_outer_phi_function(t[i])))
            A_inner_theta_list.append(om_inner_theta_function(t[i])/(om_inner_phi_function(t[i])-om_outer_phi_function(t[i])))
            A_outer_r_list.append(om_outer_r_function(t[i])/(om_inner_phi_function(t[i])-om_outer_phi_function(t[i])))
            A_outer_theta_list.append(om_outer_theta_function(t[i])/(om_inner_phi_function(t[i])-om_outer_phi_function(t[i])))
                
        
        #Now we start the big check
        label.tolist()
        [int(num) for num in label]
        for n_inner in range(int_min,int_max): 
            for n_outer in range(int_min,int_max):
                for k_inner in range(int_min,int_max):
                    for k_outer in range(int_min,int_max):

                        #Selection rule
                        if abs(k_inner - k_outer) % 2 == 0:

                            #Calculating potential m's
                            m = []
                            for j in range(len(A_inner_r_list)):
                                m.append(-n_inner*A_inner_r_list[j]-k_inner*A_inner_theta_list[j]+n_outer*A_outer_r_list[j]+k_outer*A_outer_theta_list[j])

                            #Finding where there is a resonance crossing
                            for i in range(len(m)-1):
                                if math.floor(m[i+1]) != math.floor(m[i]):
                                    if abs(m[i])<=2: # High m fourier modes will be insignificant.
                                        m_i=int(m[i])
                                        #np.gcd<=1

                                        #Screening out fundamental resonances
                                        gd = np.gcd.reduce([n_inner, n_outer, k_inner, k_outer, m_i])
                                        if gd < 3:
                                            k_checker = ((k_inner-k_outer)/2) % 2
                                            if gd == 1 or k_checker != 0:

                                                if t[i] != 1.0:
                                                    time_value_resonance = t[i]

                                                    #Final calculation and print out
                                                    for index in range(len(label)-1):
                                                        if time_value_1[index] < time_value_resonance and time_value_1[index+1] > time_value_resonance:
                                                            
                                                            #Calculates th
                                                            omega_inner_after = om_inner_r_list[index+1]*n_inner + om_inner_theta_list[index+1]*k_inner + om_inner_phi_list[index+1]*(math.floor(m[i+1]))
                                                            omega_inner_before = om_inner_r_list[index]*n_inner + om_inner_theta_list[index]*k_inner + om_inner_phi_list[index]*(math.floor(m[i+1]))
                                                            delta_t_inner = delta_t_list_inner[index+1]
                                                            omega_inner_dot = (omega_inner_after-omega_inner_before)/delta_t_inner

                                                            omega_outer_after = om_outer_r_list[index+1]*n_outer + om_outer_theta_list[index+1]*k_outer + om_outer_phi_list[index+1]*(math.floor(m[i+1]))
                                                            omega_outer_before = om_outer_r_list[index]*n_outer + om_outer_theta_list[index]*k_outer + om_outer_phi_list[index]*(math.floor(m[i+1]))
                                                            delta_t_outer = delta_t_list_outer[index+1]
                                                            omega_outer_dot = (omega_outer_after-omega_outer_before)/delta_t_outer

                                                            gamma = - mu_inner*omega_inner_dot + mu_outer*omega_outer_dot

                                                            J_inner = [J_inner_r_function(t[i]), J_inner_theta_function(t[i]), J_inner_phi_function(t[i])]
                                                            J_outer = [J_outer_r_function(t[i]), J_outer_theta_function(t[i]), J_outer_phi_function(t[i])]

                                                            _, EQL_outer = ckerr_j2eql(J_outer, M, astar)
                                                            _, EQL_inner = ckerr_j2eql(J_inner, M, astar)

                                                            _, _, anc_outer = ckerr_eql2j(list(EQL_outer), M, astar)
                                                            _, _, anc_inner = ckerr_eql2j(list(EQL_inner), M, astar)

                                                            #Calculates the values of the omegas at the exact instant of resonance
                                                            _, M_outer = ckerr_minverse(J_outer, M, astar)
                                                            omega_outer = ckerr_minv2omega(M_outer)
                                                            _, M_inner = ckerr_minverse(J_inner, M, astar)
                                                            omega_inner = ckerr_minv2omega(M_inner)

                                                            #Outputting to text file
                                                            number = number + 1
                                                            new_file_num = file_number + 1
                                                            f.write(f"{number} {new_file_num} {t[i]} {omega_inner_dot} {omega_outer_dot} {mu_inner} {mu_outer} {gamma} {m_i} {n_inner} {k_inner} {n_outer} {k_outer} {J_inner[0]} {J_inner[1]} {J_inner[2]} {J_outer[0]} {J_outer[1]} {J_outer[2]} {omega_inner[0]} {omega_inner[1]} {omega_inner[2]} {omega_outer[0]} {omega_outer[1]} {omega_outer[2]} {anc_inner[0]} {anc_inner[1]} {anc_inner[2]} {anc_outer[0]} {anc_outer[1]} {anc_outer[2]}\n")
            
                                            else:
                                                continue
            
