from cmath import sqrt
from json.tool import main
import random
import math
from gr_wrapper import ckerr_j2eql, ckerr_eql2j, rk4_j2jdot, ckerr_minverse, ckerr_minv2omega




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



#Inputs for the function
print("Initial function values:\n")
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
print("Mode vector values:\n")
print("Enter n outer:")
n_outer = float(input())
print("Enter k outer:")
k_outer = float(input())
print("Enter n inner:")
n_inner = float(input())
print("Enter k inner:")
k_inner = float(input())
print("Enter m:")
m = float(input())
print("rk4 values:\n")
print("Enter time-step:")
dt = float(input())
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

   #Running rk4
   J_r_outer_list, J_phi_outer_list, J_theta_outer_list = rk4_j2jdot(dt, t0, num_steps, J_outer[0], J_outer[1], J_outer[2])
   J_r_inner_list, J_phi_inner_list, J_theta_inner_list = rk4_j2jdot(dt, t0, num_steps, J_inner[0], J_inner[1], J_inner[2])

   #Omega calculation for each J evolution
   for j in range(num_steps):

      #Assigning J's to a list
      J_evol_outer = [J_r_outer_list[j], J_phi_outer_list[j], J_theta_outer_list[j]]
      J_evol_inner = [J_r_inner_list[j], J_phi_inner_list[j], J_theta_inner_list[j]]

      #Conversions of evolved Js to Omegas
      _, Minv_outer = ckerr_minverse(J_evol_outer, M, astar)
      _, Minv_inner = ckerr_minverse(J_evol_inner, M, astar)
      Om_outer = ckerr_minv2omega(Minv_outer)
      Om_inner = ckerr_minv2omega(Minv_inner)

      #Dot product for both sides
      dot_product_outer = n_outer * Om_outer[0] + k_outer * Om_outer[1] + m * Om_outer[2]
      dot_product_inner = n_inner * Om_inner[0] + k_inner * Om_inner[1] + m * Om_inner[2]
      print(dot_product_outer, dot_product_inner)


