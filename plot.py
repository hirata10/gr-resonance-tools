import pandas as pd
import matplotlib.pyplot as plt
import sys

# Declare if you want a plot of Omega vs J or EQL vs J
output = sys.argv[1]

def plot_columns_vs_column(file_path, J_column, Omega_r_column, Omega_theta_column, Omega_phi_column, E_column, Q_column, L_column):
    """
    Reads a text file and plots one column vs another column.

    :param file_path: Path to the text file.
    :param x_column: Column name or index for the x-axis.
    :param y_column: Column name or index for the y-axis.
    :param delimiter: Delimiter used in the text file (default is tab).
    """
    # Read the text file into a DataFrame
    df = pd.read_csv(file_path, delim_whitespace=True)

    print("Available columns in the file:")
    print(df.columns.tolist())

    # Enter the columns of the data file
    J_column = df[J_column].values
    Omega_r_column = df[Omega_r_column].values
    Omega_theta_column = df[Omega_theta_column].values
    Omega_phi_column = df[Omega_phi_column].values
    E_column = df[E_column].values
    Q_column = df[Q_column].values
    L_column = df[L_column].values
    

    # # Check if the columns exist in the DataFrame
    # if x_column not in df.columns:
    #     raise ValueError(f"Column '{x_column}' not found in the file.")
    # if y_column not in df.columns:
    #     raise ValueError(f"Column '{y_column}' not found in the file.")

    # Plotting
    #Have to adjust axes and title of plot manually to the appropriate J that is varying

    if output == "Omega":
        plt.figure(figsize=(10, 6))
        plt.plot(J_column, Omega_r_column, color='b', label = r"$\Omega^r$")
        plt.plot(J_column, Omega_theta_column, color='r', label = r"$\Omega^{\theta}$")
        plt.plot(J_column, Omega_phi_column, color='g', label = r"$\Omega^{\phi}$")
        #plt.plot(df[x_column], df[y_column], marker='o', linestyle='-', color='b')
        #plt.scatter(df[x_column], df[y_column], color='blue', marker='o')
        plt.xlabel(r"$J_{\theta}$")
        plt.ylabel(r"$\Omega^i$")
        plt.title(r'Orbital frequencies vs $J_{\theta}$')
        plt.legend(loc="lower left")
        plt.grid(True)
        plt.show()
    
    elif output == "EQL":
        plt.figure(figsize=(10, 6))
        plt.plot(J_column, E_column, color='b', label = r"E")
        plt.plot(J_column, Q_column, color='r', label = r"Q")
        plt.plot(J_column, L_column, color='g', label = r"L")
        #plt.plot(df[x_column], df[y_column], marker='o', linestyle='-', color='b')
        #plt.scatter(df[x_column], df[y_column], color='blue', marker='o')
        plt.xlabel(r"$J_{\theta}$")
        plt.ylabel(r"EQL")
        plt.title(r'EQL vs $J_{\theta}$')
        plt.legend(loc="lower left")
        plt.grid(True)
        plt.show()
    
    else:
        print("Did not specify proper output option: Omega or EQL \n")

# Example usage
file_path = 'J2Omega_vary_phi.txt'  # Replace with your file path
J_column = "J_phi"  # Replace with the column for the x-axis
Omega_r_column = "Omega_r"  # Replace with the column for the y-axis
Omega_theta_column = "Omega_theta"
Omega_phi_column = "Omega_phi"
E_column = "E"  # Replace with the column for the y-axis
Q_column = "Q"
L_column = "L"

plot_columns_vs_column(file_path, J_column, Omega_r_column, Omega_theta_column, Omega_phi_column, E_column, Q_column, L_column)
