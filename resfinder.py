
from re import M
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
from sympy import Point, Line
from scipy import interpolate
import sys
from gr_wrapper import ckerr_minverse, ckerr_minv2omega, ckerr_eql2j, ckerr_j2eql
import os             
import math




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



#Note that in the paper, we chose 0.00001, 0.00001, 1, 0.9, and ran over 50 files
print("Enter the EMRI body's mass:")
mu_inner = 0.00001 #float(input())
print("Enter the perturber's mass:")
mu_outer = 0.00001 #float(input())
print("Enter the SMBH mass:")
M = 1 #float(input())
print("Enter the spin parameter:")
astar = 0.9 #float(input())
print("Enter the number of inner/outer body file pairs:")
file_amount = 50 #float(input())



#Starting the runs
number = 0
for file_number in range(file_amount):

    #Importing all the info from the files. Feel free to change this depending on the files you made
    #Also note how we use label, because there will (very likely) always be less points in outer than inner.
    _, time_value_1, J_r_inner_list, J_theta_inner_list, J_phi_inner_list, om_inner_r_list, om_inner_theta_list, om_inner_phi_list, delta_t_list_inner = np.loadtxt("inner_body_runs/J_evolve_inner_" + str(file_number+1) + ".txt", unpack=True)
    label, time_value_2, J_r_outer_list, J_theta_outer_list, J_phi_outer_list, om_outer_r_list, om_outer_theta_list, om_outer_phi_list, delta_t_list_outer = np.loadtxt("outer_body_runs/J_evolv_outer_" + str(file_number+1) + ".txt", unpack=True) 
    
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
    #Note (for the paper and in general) we chose to range between -5 and 5 bc anything greater, we found, was insignificant 
    for n_inner in range(-5,5): 
        for n_outer in range(-5,5):
            for k_inner in range(-5,5):
                for k_outer in range(-5,5):

                    #Selection rule
                    if abs(k_inner - k_outer) == 2:
                        #Calculating potential m's
                        m = []
                        for j in range(len(A_inner_r_list)):
                            m.append(-n_inner*A_inner_r_list[j]-k_inner*A_inner_theta_list[j]+n_outer*A_outer_r_list[j]+k_outer*A_outer_theta_list[j])

                        #Finding where there is a resonance crossing
                        for i in range(len(m)-1):
                            if math.floor(m[i+1]) != math.floor(m[i]):
                                if abs(m[i])<=2:
                                    if t[i] != 1.0:
                                        time_value_resonance = t[i]

                                        #Final calculation and print out
                                        for index in range(len(label)-1):
                                            if time_value_1[index] < time_value_resonance and time_value_1[index+1] > time_value_resonance:
                                                
                                                omega_inner_after = om_inner_r_list[index+1]*n_inner + om_inner_theta_list[index+1]*k_inner + om_inner_phi_list[index+1]*(math.floor(m[i+1]))
                                                omega_inner_before = om_inner_r_list[index]*n_inner + om_inner_theta_list[index]*k_inner + om_inner_phi_list[index]*(math.floor(m[i+1]))
                                                delta_t_inner = delta_t_list_inner[index+1]
                                                omega_inner_dot = (omega_inner_after-omega_inner_before)/delta_t_inner

                                                omega_outer_after = om_outer_r_list[index+1]*n_outer + om_outer_theta_list[index+1]*k_outer + om_outer_phi_list[index+1]*(math.floor(m[i+1]))
                                                omega_outer_before = om_outer_r_list[index]*n_outer + om_outer_theta_list[index]*k_outer + om_outer_phi_list[index]*(math.floor(m[i+1]))
                                                delta_t_outer = delta_t_list_outer[index+1]
                                                omega_outer_dot = (omega_outer_after-omega_outer_before)/delta_t_outer

                                                gamma = mu_inner*omega_inner_dot-mu_outer*omega_outer_dot

                                                J_inner = [J_inner_r_function(t[i]), J_inner_theta_function(t[i]), J_inner_phi_function(t[i])]
                                                J_outer = [J_outer_r_function(t[i]), J_outer_theta_function(t[i]), J_outer_phi_function(t[i])]

                                                _, EQL_outer = ckerr_j2eql(J_outer, M, astar)
                                                _, EQL_inner = ckerr_j2eql(J_inner, M, astar)

                                                _, _, anc_outer = ckerr_eql2j(list(EQL_outer), M, astar)
                                                _, _, anc_inner = ckerr_eql2j(list(EQL_inner), M, astar)

                                                number = number + 1

                                                print(number, file_number + 1, t[i], omega_inner_dot, omega_outer_dot, mu_inner, mu_outer, gamma, math.floor(m[i+1]), n_inner, k_inner, n_outer, k_outer, J_inner[0], J_inner[1], J_inner[2], J_outer[0], J_outer[1], J_outer[2], anc_inner[0], anc_inner[1], anc_inner[2], anc_outer[0], anc_outer[1], anc_outer[2])
    
 

    
