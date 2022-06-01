from cmath import sqrt
from json.tool import main
import random
import math
from gr_wrapper import ckerr_j2eql, ckerr_eql2j, ckerr_minverse, ckerr_minv2omega




"""
Hello! Here, we have three random number
generators set to each of the action variables 
as well as a series of checks to get viable
actions for future steps. Note these checks
are valid in a Keplerian model (scaling GM=1). 

Inputs: outer/inner sm axis, spin parameter, mass 
of central body, number of resonance pairs needed.

Outputs: List of Omega vector pairs for subsequent
rk4 integration/other purposes.
"""



#Inputs for the function
print("Enter the outer semi-major axis:")
sm_axis_outer = float(input())
print("Enter the inner semi-major axis:")
sm_axis_inner = float(input())
print("Enter spin parameter (from 0 to 1):")
astar = float(input()) 
print("Enter mass of the central body:")
M = float(input())
print("Enter number of resonance pairs needed:")
max = float(input())
print("Calculation starting. Grab some popcorn!")


#Variable assignments
sqrt_sm_axis_outer = math.sqrt(sm_axis_outer)
sqrt_sm_axis_inner = math.sqrt(sm_axis_inner)

rk4_list_outer = []
rk4_list_inner = []




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

   #Check 1: sum(J) = sqrt of sm-axes
   if J_r_outer + J_theta_outer + J_phi_outer != sqrt_sm_axis_outer or J_r_inner + J_theta_inner + J_phi_inner != sqrt_sm_axis_inner:
      continue

   #Check 2: Closed orbit (EQL[0] < 1)
   _, EQL_outer = ckerr_j2eql(J_outer, M, astar)
   _, EQL_inner = ckerr_j2eql(J_inner, M, astar)
   if EQL_outer[0] > 1 or EQL_inner[0] > 1:
      continue

   #Check 3: r_peri_outer > 3 r_apo_inner
   anc_outer = ckerr_eql2j(EQL_outer, M, astar)
   anc_inner = ckerr_eql2j(EQL_inner, M, astar)
   if anc_outer[2] < 3*anc_inner[1]:
      continue

   #Conversion to omegas
   _, Minv_outer = ckerr_minverse(J_outer, M, astar)
   _, Minv_inner = ckerr_minverse(J_inner, M, astar)
   Om_outer = ckerr_minv2omega(Minv_outer)
   Om_inner = ckerr_minv2omega(Minv_inner)


   rk4_list_outer += [J_outer]
   rk4_list_inner += [J_inner]
   i += 1

print(rk4_list_outer, rk4_list_inner)

