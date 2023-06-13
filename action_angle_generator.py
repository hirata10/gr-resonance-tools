from cmath import sqrt
from json.tool import main
import time
import random
import math
from gr_wrapper import ckerr_j2eql, ckerr_eql2j, rk4_j2jdot, ckerr_minverse, ckerr_minv2omega, j2jdot_component



"""
Hello! Here, we have three random number
generators set to each of the action variables 
as well as a series of checks to get viable
actions for further analysis. Note these checks
are valid in a Keplerian model (scaling GM=1). 

Inputs: outer/inner sm axis, spin parameter, mass 
of central body, number of resonance pairs needed.

Outputs: Inner J's, outer J's, and ancilliary arrays.
"""



print("Enter the outer semi-major axis:")
sm_axis_outer = 300 #float(input())
print("Enter the inner semi-major axis:")
sm_axis_inner = 50 #float(input())
print("Enter mass of the central body:")
M = 1.0 #float(input())
print("Enter the spin parameter of the central body:")
astar = 0.9 #float(input())
print("Enter mass of outer body:")
mu_body_outer = 1 #float(input())
print("Enter mass of inner body:")
mu_body_inner = 1 #float(input())
print("Enter number of resonance pairs needed:")
max = 50 #int(input())
print("Enter initial time:")
t0 = 1 #float(input())
print("Enter number of steps:")
num_steps = 20 #int(input())
print("Calculation starting. Grab some popcorn!")


#Variable assignments
sqrt_sm_axis_outer = math.sqrt(sm_axis_outer)
sqrt_sm_axis_inner = math.sqrt(sm_axis_inner)
r_inner_solution = M*(1+math.sqrt(1-astar**2))


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
   #print("J_inner", J_r_inner, J_theta_inner, J_phi_inner)
   #print("J_outer", J_r_outer, J_theta_outer, J_phi_outer)
   #print(J_r_inner, J_theta_inner, J_phi_inner, anc_inner, J_r_outer, J_theta_outer, J_phi_outer, anc_outer, "Starting checks!")


   #Check 1: sum(J) = sqrt of sm-axes
   if J_r_outer + J_theta_outer + J_phi_outer != sqrt_sm_axis_outer or J_r_inner + J_theta_inner + J_phi_inner != sqrt_sm_axis_inner:
      #print("Check 1 failed!")
      continue

   #Check 2: Closed orbit (EQL[0] < 1)
   outer_checker, EQL_outer = ckerr_j2eql(J_outer, M, astar)
   inner_checker, EQL_inner = ckerr_j2eql(J_inner, M, astar)
   if EQL_outer[0] > 1 or EQL_inner[0] > 1:
      #print("Check 2a failed!")
      continue
   #if outer_checker != 1 or inner_checker != 1: 
      #print("Check 2b failed!")
      #continue

   #Check 3: r_peri_outer > 3 r_apo_inner
   _, _, anc_outer = ckerr_eql2j(list(EQL_outer), M, astar)
   _, _, anc_inner = ckerr_eql2j(list(EQL_inner), M, astar)
   if anc_outer[1] < 2*anc_inner[2]:
      #print("Check 3 failed!")
      continue

   #Check 4: r_peri_outer, r_peri_inner > 2M
   if anc_outer[1] < r_inner_solution or anc_inner[1] < r_inner_solution:
      #print("Check 4 failed!")
      continue
   
   
   #Print statement and increment
   i += 1
   print(i, J_r_inner, J_theta_inner, J_phi_inner, J_r_outer, J_theta_outer, J_phi_outer, anc_inner, anc_outer)
   