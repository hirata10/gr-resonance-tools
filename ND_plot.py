import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

output = sys.argv[1]

def compute_derivative(file_path, J_column, Omega_r_column, Omega_theta_column, Omega_phi_column, E_column, Q_column, L_column):
    # Read the data from the text file
    data = pd.read_csv(file_path, delim_whitespace=True)
    print("Available columns in the file:")
    print(data.columns.tolist())

    # # Check if the specified columns exist
    # if column_x not in data.columns or column_y not in data.columns:
    #     raise ValueError(f"Columns '{column_x}' or '{column_y}' not found in the data file.")

    # Extract the data for the specified columns
    J = data[J_column].values
    Omega_r = data[Omega_r_column].values
    Omega_theta = data[Omega_theta_column].values
    Omega_phi = data[Omega_phi_column].values
    E = data[E_column].values
    Q = data[Q_column].values
    L = data[L_column].values

    # Compute the numerical derivative
    dOmega_r_dJ = np.gradient(J, Omega_r)
    dOmega_theta_dJ = np.gradient(J, Omega_theta)
    dOmega_phi_dJ = np.gradient(J, Omega_phi)
    dE_dJ = np.gradient(J, E)
    dQ_dJ = np.gradient(J, Q)
    dL_dJ = np.gradient(J, L)

    return J, dOmega_r_dJ, dOmega_theta_dJ, dOmega_phi_dJ, dE_dJ, dQ_dJ, dL_dJ

#Have to adjust axes and title of plot manually to the appropriate J that is varying
def plot_derivative(J, dOmega_r_dJ, dOmega_theta_dJ, dOmega_phi_dJ, output):
    if output == "Omega":
        plt.figure(figsize=(10, 6))
        plt.plot(J, dOmega_r_dJ, label=r'd$\Omega^r$/d$J_{\theta}$', color='blue')
        plt.plot(J, dOmega_theta_dJ, label=r'd$\Omega^{\theta}$/d$J_{\theta}$', color='red')
        plt.plot(J, dOmega_phi_dJ, label=r'd$\Omega^{\phi}$/d$J_{\theta}$', color='green')
        plt.xlabel(r"$J_{\theta}$")
        plt.ylabel(r'd$\Omega^i$/d$J_{\theta}$')
        plt.title(r'Derivative of orbital frequencies with respect to $J_{\theta}$')
        plt.legend()
        plt.grid()
        plt.show()
    
    elif output == "EQL":
        plt.figure(figsize=(10, 6))
        plt.plot(J, dE_dJ, label=r'dE/d$J_{\theta}$', color='blue')
        plt.plot(J, dQ_dJ, label=r'dQ/d$J_{\theta}$', color='red')
        plt.plot(J, dL_dJ, label=r'dL/d$J_{\theta}$', color='green')
        plt.xlabel(r"$J_{\theta}$")
        plt.ylabel(r'dEQL/d$J_{\theta}$')
        plt.title(r'Derivative of EQL with respect to $J_{\theta}$')
        plt.legend()
        plt.grid()
        plt.show()
    
    else:
        print("Did not specify valid output: Omega or EQL \n")


# Specify the file path and columns for differentiation
file_path = 'J2Omega_vary_phi.txt'  # Change this to your data file with the appropriate changing J
J_column = 'J_phi'    # Change this to whichever J component is changing (the others should be fixed and it should match the data file used above)
Omega_r_column = 'Omega_r' 
Omega_theta_column = 'Omega_theta'
Omega_phi_column = 'Omega_phi'
E_column = 'E'
Q_column = 'Q'
L_column = 'L'

J, dOmega_r_dJ, dOmega_theta_dJ, dOmega_phi_dJ, dE_dJ, dQ_dJ, dL_dJ = compute_derivative(file_path, J_column, Omega_r_column, Omega_theta_column, Omega_phi_column, E_column, Q_column, L_column)
#print(J, dL_dJ)
plot_derivative(J, dOmega_r_dJ, dOmega_theta_dJ, dOmega_phi_dJ, output)
