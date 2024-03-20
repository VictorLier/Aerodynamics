import os
import subprocess
import numpy as np

# %% Inputs

airfoil_name = "NACA4415"
alpha_i = -6
alpha_f = 10
alpha_step = 0.5
Re = 6e6
n_iter = 100

# %% XFOIL input file writer

if os.path.exists("polar_file.txt"):
    os.remove("polar_file.txt")

input_file = open("input_file.in", 'w')
input_file.write("LOAD {0}.dat\n".format(airfoil_name))
input_file.write(airfoil_name + '\n')
input_file.write("PANE\n")
input_file.write("OPER\n")
input_file.write("Visc {0}\n".format(Re))
input_file.write("PACC\n")
input_file.write("polar_file.txt\n\n")
input_file.write("ITER {0}\n".format(n_iter))
input_file.write("ASeq {0} {1} {2}\n".format(alpha_i, alpha_f,
                                             alpha_step))
input_file.write("\n\n")
input_file.write("quit\n")
input_file.close()

subprocess.call("xfoil.exe < input_file.in", shell=True)

polar_data = np.loadtxt("polar_file.txt", skiprows=12)

import numpy as np
import matplotlib.pyplot as plt

# Read data from the text file
file_path = 'polar_file.txt'  # Update with the actual path to your text file
data = np.genfromtxt(file_path, skip_header=11)  # Skip the first 4 lines which contain headers

# Extract alpha and CL data
alpha = data[:, 0]
CL = data[:, 1]
CD = data[:, 2]

# Plot alpha vs. CL
plt.subplot(2,1,1)
plt.plot(alpha, CL)
plt.ylabel('Coefficient of Lift (CL)')
plt.title('Angle of Attack vs. Coefficient of Lift/Drag')
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(alpha, CD)
plt.xlabel('Angle of Attack (alpha)')
plt.ylabel('Coefficient of Drag (CD)')
plt.grid(True)
plt.tight_layout()




plt.show()