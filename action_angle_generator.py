from cmath import sqrt
from json.tool import main
import time
import random
import math
from gr_wrapper import ckerr_j2eql, ckerr_eql2j
import globalpars


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

#Defining variables for random number generator
M = globalpars.GLOBALPAR_M
astar = globalpars.GLOBALPAR_astar
sm_axis_outer = globalpars.GLOBALPAR_sm_axis_outer
sm_axis_inner = globalpars.GLOBALPAR_sm_axis_inner
max = globalpars.GLOBALPAR_N_sys
sqrt_sm_axis_outer = math.sqrt(sm_axis_outer)
sqrt_sm_axis_inner = math.sqrt(sm_axis_inner)
r_inner_solution = M*(1+math.sqrt(1-astar**2))
print("SMBH mass and spin: ", M, astar)

# Open a text file for writing
with open('action_angle_pairs.txt', 'w') as f:

    # Loop to generate random actions and perform checks
    i = 0
    while i < max:
       # Generate random actions (transformation method)
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

       # Perform checks
       # Check 1: sum(J) = sqrt of sm-axes
       if J_r_outer + J_theta_outer + J_phi_outer != sqrt_sm_axis_outer or J_r_inner + J_theta_inner + J_phi_inner != sqrt_sm_axis_inner:
          continue

       # Check 2: Closed orbit (EQL[0] < 1; EQL[0] = E, EQL[1] = Q, EQL[2] = L)
       outer_checker, EQL_outer = ckerr_j2eql(J_outer, M, astar)
       inner_checker, EQL_inner = ckerr_j2eql(J_inner, M, astar)
       if EQL_outer[0] > 1 or EQL_inner[0] > 1:
          continue

       # Check 3: r_peri_outer > 2 r_apo_inner (avoid crossings); anc[0] = inclination, anc[1] = pericenter, anc[2] = apocenter
       _, _, anc_outer = ckerr_eql2j(list(EQL_outer), M, astar)
       _, _, anc_inner = ckerr_eql2j(list(EQL_inner), M, astar)
       if anc_outer[1] < 2*anc_inner[2]:
          continue

       # Check 4: r_peri_outer, r_peri_inner > 2M (outside the outer horizon of Kerr BH)
       if anc_outer[1] < r_inner_solution or anc_inner[1] < r_inner_solution:
          continue
       
       #Check 5: eccentricity < 0.8 (require greater numerical resolution in r-direction, future work)
       e_inner = (anc_inner[2] - anc_inner[1]) / (anc_inner[2] + anc_inner[1])
       e_outer = (anc_outer[2] - anc_outer[1]) / (anc_outer[2] + anc_outer[1])
       if e_inner > 0.8 or e_outer > 0.8:
          continue

       # Write output to file
       i += 1
       f.write(f"{i} {J_r_inner} {J_theta_inner} {J_phi_inner} {J_r_outer} {J_theta_outer} {J_phi_outer}\n")


#Maybe an extra check that accounts for the mass ratio or? Idk.
#Why is this not showing up in output.txt