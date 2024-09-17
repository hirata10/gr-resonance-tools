import numpy as np
import scipy.optimize as opt
#from sympy import Point, Line
from scipy import interpolate
import sys
from gr_wrapper import ckerr_eql2j, ckerr_j2eql, ckerr_minverse, ckerr_minv2omega
import os             
import math
import sys
import globalpars

M = globalpars.GLOBALPAR_M
astar = globalpars.GLOBALPAR_astar
mu_outer = globalpars.GLOBALPAR_mu_outer
mu_inner = globalpars.GLOBALPAR_mu_inner
file_amount = globalpars.GLOBALPAR_N_sys
int_min = globalpars.GLOBALPAR_N_min_Fourier
int_max = globalpars.GLOBALPAR_N_max_Fourier
N_step = 50
vary_label = int(sys.argv[1])

print(vary_label)

#Input array for action variables. Single value example
#J = [0.039287, 1.559630, 5.472151]

#Class that defines the output of a function to the txt file and 
class DualOutput:
    def __init__(self, file_name):
        self.file = open(file_name, 'w')  # Open the file for writing
        self.terminal = sys.stdout  # Save a reference to the original stdout
    
    def write(self, message):
        self.file.write(message)  # Write to the file
        self.terminal.write(message)  # Write to the terminal
    
    def flush(self):
        self.file.flush()  # Ensure the file buffer is flushed
        self.terminal.flush()  # Ensure the terminal buffer is flushed

    def close(self):
        self.file.close()  # Close the file when done

# filename = "J2Omega/J2Omega.txt"

# directory = os.path.dirname(filename)
# if not os.path.exists(directory) and directory != '':
#     os.makedirs(directory)

#Function that maps the action variables to the orbital frequencies and EQL using functions from kerrtraj.c
def pyJ2Omega(input_J, mass, spin):
    _, EQL = ckerr_j2eql(input_J, mass, spin)

    _, M = ckerr_minverse(input_J, mass, spin)

    omega = ckerr_minv2omega(M)

    return(omega, EQL)


#dual_output = DualOutput('J2Omega_vary_phi.txt')
    
# Redirect sys.stdout to the DualOutput instance
#sys.stdout = dual_output
#print("J_phi vs Omega_i")
#print("J_r \t J_theta \t J_phi \t Omega_r \t Omega_theta \t Omega_phi \t E \t Q \t L \n")
### For loop for each component of J changes while keeping the others fixed
# for i in range(N_step): 
#     for j in range(N_step):
#         for k in range(N_step):
#             J = [0.04 + i/N_step, 1.55 + j/N_step, 5.47 + k/N_step]
#             #J = [0.039287, 1.559630, 5.472151]
#             Omega, EQL = pyJ2Omega(J, M, astar)
#             print(J[0], J[1], J[2], Omega[0], Omega[1], Omega[2], EQL[0], EQL[1], EQL[2], end='\n')  # Print the output, end='' prevents double new lines

### For loops where a single component of J is varied while the other two remain fixed. Ideal for computing numerical derivatives of Omega^j w.r.t J_i

### Vary J_phi
if (vary_label == 2):
    print("please?")

    dual_output = DualOutput('J2Omega_vary_phi.txt')
    # Redirect sys.stdout to the DualOutput instance
    sys.stdout = dual_output
    print("J_phi vs Omega_i")
    print("J_r \t J_theta \t J_phi \t Omega_r \t Omega_theta \t Omega_phi \t E \t Q \t L \n")
    for i in range(N_step):
        J = [0.04, 1.55, 5.47 + i/N_step]
        #J = [0.039287, 1.559630, 5.472151]
        Omega, EQL = pyJ2Omega(J, M, astar)
        print(J[0], J[1], J[2], Omega[0], Omega[1], Omega[2], EQL[0], EQL[1], EQL[2], end='\n')  # Print the output, end='' prevents double new lines
    
     # Restore original stdout
    sys.stdout = dual_output.terminal

    # Close the file when done
    dual_output.close()
    

### Vary J_theta
elif (vary_label == 1):
    dual_output = DualOutput('J2Omega_vary_theta.txt')
    
    # Redirect sys.stdout to the DualOutput instance
    sys.stdout = dual_output
    print("J_theta vs Omega_i")
    print("J_r \t J_theta \t J_phi \t Omega_r \t Omega_theta \t Omega_phi \t E \t Q \t L \n")
    for i in range(N_step):
        J = [0.04, 1.55 + i/N_step, 5.47]
        #J = [0.039287, 1.559630, 5.472151]
        Omega, EQL = pyJ2Omega(J, M, astar)
        print(J[0], J[1], J[2], Omega[0], Omega[1], Omega[2], EQL[0], EQL[1], EQL[2], end='\n')  # Print the output, end='' prevents double new lines
    # Restore original stdout
    sys.stdout = dual_output.terminal

    # Close the file when done
    dual_output.close()


### Vary J_r
elif (vary_label == 0):
    dual_output = DualOutput('J2Omega_vary_r.txt')
    
    # Redirect sys.stdout to the DualOutput instance
    sys.stdout = dual_output
    print("J_r vs Omega_i")
    print("J_r \t J_theta \t J_phi \t Omega_r \t Omega_theta \t Omega_phi \t E \t Q \t L \n")
    for i in range(N_step):
        J = [0.04 + i/N_step, 1.55, 5.47]
        #J = [0.039287, 1.559630, 5.472151]
        Omega, EQL = pyJ2Omega(J, M, astar)
        print(J[0], J[1], J[2], Omega[0], Omega[1], Omega[2], EQL[0], EQL[1], EQL[2], end='\n')  # Print the output, end='' prevents double new lines

    # Restore original stdout
    sys.stdout = dual_output.terminal

    # Close the file when done
    dual_output.close()

else:
    print("Invalid J component label")

# # Restore original stdout
# sys.stdout = dual_output.terminal

# # Close the file when done
# dual_output.close()

#Omega= pyJ2Omega(J, M, astar)



# print("Orbtial frequencies:", omega)
# print("Adiabatic constants of motion:", EQL)