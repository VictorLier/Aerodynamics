import os
import subprocess
from numpy import pi, sin, cos, exp, sqrt, log
import matplotlib.pyplot as plt
import numpy as np
airfoil_name = "NACA6412"
Re = 1000000
n_iter = 100

# Xfoil parameters
alpha = 10
Re = 1000000
n_iter = 300

# %% XFOIL input file writer
input_file = open("input_file.in", 'w')
input_file.write("LOAD {0}.dat\n".format(airfoil_name))
input_file.write(airfoil_name + '\n')
input_file.write("OPER\n")
input_file.write("Visc {0}\n".format(Re))
input_file.write("ITER {0}\n".format(n_iter))
input_file.write("ALFA {}\n".format(alpha))  # Add newline character here
input_file.write("CPWR cp_out.txt\n")  # Add newline character here
input_file.write("QUIT\n")
input_file.close()

# Run Xfoil
subprocess.call("xfoil.exe < input_file.in", shell=True)

# Check if Xfoil has generated the output file
CP_data = np.loadtxt("cp_out.txt", skiprows=4)

x = CP_data[:,0]
CP = CP_data[:,2]

plt.plot(x, -CP)
plt.xlabel("Chord")
plt.ylabel("- CP")
plt.title("- Coefficient of pressure")
plt.suptitle("$\\alpha = {}$".format(alpha))
plt.grid()
plt.show()
