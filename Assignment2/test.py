import numpy as np

nA = 10
dpi = np.pi/nA
theta_array = np.linspace(dpi, np.pi-dpi, nA)

alpha_L0 = -0.0854964326076398
alpha = 0.17453292519943295

m0 = 5.856774581827276
b = 7853.981633974483

AR = 4

# Initialize ordinate vector - o
o = np.array([])

# Initialize coefficient matrix - B
B = np.zeros((nA, nA))

# Fill in B matrix
for j, theta in enumerate(theta_array):
    for i in range(nA):
        B[j, i] = (-4*AR/m0 * np.sin((1+i)*theta) - (1+i)*(np.sin((1+i)*theta)/np.sin(theta))) # Silde 26 week 8
    o = np.append(o, alpha - alpha_L0)
# Solve for A
A = np.linalg.solve(B, o)

print("A", np.round(A,7))