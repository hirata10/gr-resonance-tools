from cmath import sqrt
from json.tool import main
import random
import math
from gr_wrapper import ckerr_j2eql, ckerr_eql2j, rk4_j2jdot, ckerr_minverse, ckerr_minv2omega, j2jdot_component




"""
Hello! Here, we have three random number
generators set to each of the action variables 
as well as a series of checks to get viable
actions for future steps. Note these checks
are valid in a Keplerian model (scaling GM=1). 

Inputs: outer/inner sm axis, spin parameter, mass 
of central body, number of resonance pairs needed.

Outputs: Dot products of evolved Omega's and mode
vectors after running through rk4.
"""



print("Enter the outer semi-major axis:")
sm_axis_outer = float(input())
print("Enter the inner semi-major axis:")
sm_axis_inner = float(input())
print("Enter mass of the central body:")
M = float(input())
print("Enter mass of outer body:")
mu_body_outer = float(input())
print("Enter mass of inner body:")
mu_body_inner = float(input())
print("Enter number of resonance pairs needed:")
max = int(input())
print("Enter initial time:")
t0 = float(input())
print("Enter number of steps:")
num_steps = int(input())
print("Calculation starting. Grab some popcorn!")


#Variable assignments
sqrt_sm_axis_outer = math.sqrt(sm_axis_outer)
sqrt_sm_axis_inner = math.sqrt(sm_axis_inner)


i = 0
while i < max:
   #Generate random actions
   zeta1 = random.random()
   eta1 = random.random()
   zeta2 = random.random()
   eta2 = random.random()
   u1 = math.sqrt(zeta1) * eta1
   u2 = math.sqrt(zeta1) * (1 - eta1)
   u3 = math.sqrt(zeta2) * eta2
   u4 = math.sqrt(zeta2) * (1 - eta2)
   J_r_outer = u1 * sqrt_sm_axis_outer 
   J_theta_outer = u2 * sqrt_sm_axis_outer 
   J_phi_outer = (1 - u1 - u2) * sqrt_sm_axis_outer 
   J_r_inner = u3 * sqrt_sm_axis_inner 
   J_theta_inner = u4 * sqrt_sm_axis_inner 
   J_phi_inner = (1 - u3 - u4) * sqrt_sm_axis_inner 
   J_outer = [J_r_outer, J_theta_outer, J_phi_outer]
   J_inner = [J_r_inner, J_theta_inner, J_phi_inner]
   print("J_inner", J_r_inner, J_theta_inner, J_phi_inner)
   print("J_outer", J_r_outer, J_theta_outer, J_phi_outer)


   #Check 1: sum(J) = sqrt of sm-axes
   if J_r_outer + J_theta_outer + J_phi_outer != sqrt_sm_axis_outer or J_r_inner + J_theta_inner + J_phi_inner != sqrt_sm_axis_inner:
      print("Check 1 failed!")
      continue
   

   #Check 2: Closed orbit (EQL[0] < 1)
   a, EQL_outer_01 = ckerr_j2eql(J_outer, M, 0.1)
   b, EQL_inner_01 = ckerr_j2eql(J_inner, M, 0.1)
   c, EQL_outer_09 = ckerr_j2eql(J_outer, M, 0.9)
   d, EQL_inner_09 = ckerr_j2eql(J_inner, M, 0.9)
   if EQL_outer_01[0] > 1 or EQL_inner_01[0] > 1:
      if EQL_outer_09[0] > 1 or EQL_inner_09[0] > 1:
         print("Check 2a failed!")
         continue
   if a != 1 or b != 1: 
      print("Check 2b failed!")
      continue
   if c != 1 or d != 1: 
      print("Check 2c failed!")
      continue

   #Check 3: r_peri_outer > 3 r_apo_inner
   anc_outer_01 = ckerr_eql2j(list(EQL_outer_01), M, 0.1)
   anc_inner_01 = ckerr_eql2j(list(EQL_outer_01), M, 0.1)
   anc_outer_09 = ckerr_eql2j(list(EQL_outer_09), M, 0.9)
   anc_inner_09 = ckerr_eql2j(list(EQL_outer_09), M, 0.9)
   if anc_outer_01[2] < 3*anc_inner_01[1]:
      if anc_outer_09[2] < 3*anc_inner_09[1]:
         print("Check 3 failed!")
         continue

   print("All checks passed!")
   print("Checked J_inner", J_r_inner, J_theta_inner, J_phi_inner)
   print("Checked J_outer", J_r_outer, J_theta_outer, J_phi_outer)
   _, J_dot_r, J_dot_theta, J_dot_phi = j2jdot_component(1, 1, 1, 1, J_r_inner, J_theta_inner, J_phi_inner, M, 0.1)
   timescale = J_phi_inner / J_dot_phi[0]
   print(timescale)
   print("Enter timesteps:\n")

   dt_01 = float(input())
   dt_09 = float(input())

   #Running rk4
   J_r_inner_list, J_theta_inner_list, J_phi_inner_list = rk4_j2jdot(dt_01, t0, num_steps, J_inner[0], J_inner[1], J_inner[2], mu_body_inner, M, 0.1)
   print(J_r_inner_list, J_theta_inner_list, J_phi_inner_list, 0.1, "first run done")
   J_r_outer_list, J_theta_outer_list, J_phi_outer_list = rk4_j2jdot(dt_01, t0, num_steps, J_outer[0], J_outer[1], J_outer[2], mu_body_outer, M, 0.1)
   print(J_r_outer_list, J_theta_outer_list, J_phi_outer_list, 0.9, "second run done")
   J_r_inner_list, J_theta_inner_list, J_phi_inner_list = rk4_j2jdot(dt_09, t0, num_steps, J_inner[0], J_inner[1], J_inner[2], mu_body_inner, M, 0.9)
   print(J_r_inner_list, J_theta_inner_list, J_phi_inner_list, 0.1, "third run done")
   J_r_outer_list, J_theta_outer_list, J_phi_outer_list = rk4_j2jdot(dt_09, t0, num_steps, J_outer[0], J_outer[1], J_outer[2], mu_body_outer, M, 0.9)
   print(J_r_inner_list, J_theta_inner_list, J_phi_inner_list, 0.9, "fourth run done")


   i += 1












