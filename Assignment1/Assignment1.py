import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# Parameters defining the shape of the airfoil
c = 1.2
s1, s2 = c*0.1, c*0.1
U0 = 1

n_points = 101
theta = np.linspace(0, 2*np.pi, n_points)

R = lambda theta, s1, s2: s2*np.sin(theta) - s1*np.cos(theta) + np.sqrt((s1*np.cos(theta)-s2*np.sin(theta))**2 + c*(c+2*s1))
zeta = lambda z: z + c**2/z


# Symbolic f(z)
z_s, s_s = sp.symbols('z, s', complex=True)
a_s, c_s, alpha_s, U_s, s1_s, s2_s, theta_s, Gamma_s = sp.symbols('a, c, alpha, U_\infty, s_1, s_2, theta, Gamma', real=True, positive=True)
R_s = s2_s*sp.sin(theta_s) - s1_s*sp.cos(theta_s) + sp.sqrt((s1_s*sp.cos(theta_s)-s2_s*sp.sin(theta_s))**2 + c_s*(c_s+2*s1_s))
#s1_s, s2_s = sp.re(s_s), sp.im(s_s)
Gamma_s = 4*sp.pi*U_s*a_s*sp.sin(alpha_s + sp.atan(sp.im(s_s)/(c_s-sp.re(s_s))))
f_s = U_s*( (z_s-s_s)*sp.exp(-sp.I*alpha_s) + a_s**2/(z_s-s_s)*sp.exp(sp.I*alpha_s) ) + sp.I*Gamma_s/(2*sp.pi)*sp.log((z_s-s_s)/a_s)

print("F(z) =", f_s)

f = sp.lambdify((z_s, s_s, alpha_s, a_s), f_s.subs({c_s:c, U_s:U0}), "numpy")

# Plot the original circle and airfoil
z_c = np.exp(1j*theta)
z_a = np.array([R(t,s1,s2) for t in theta])*np.exp(1j*theta)

# Find the lenght of the airfoil, the maximum thickness and the maximum camber
l = np.real( zeta(z_a)[0] - zeta(z_a)[n_points//2])
t = np.max( [np.imag(zeta(z_a)[i] - zeta(z_a)[n_points-1-i]) for i in range(n_points//2)])

h_x = np.array( [np.real(zeta(z_a)[i]) for i in range(n_points//2)])
h_y = np.array( [np.imag(zeta(z_a)[i] + zeta(z_a)[n_points-1-i])/2 for i in range(n_points//2)])

h = np.max(h_y)

print(f"s1/c = {s1/c}, s2/c = {s2/c}")
print(f"l = {l}, t = {t}, h = {h}")
print(f"t/l = {round(t/l*100,2)}, h/l = {round(h/l*100,2)}")

sym_s1 = np.array([0.0836, 0.1833, 0.3057])
sym_s2 = np.array([0, 0, 0])

asym_s1 = np.array([0.083,0.181,0.301])
asym_s2 = np.array([0.216,0.260,0.261])

Titles = ["Symmetric airfoils", "Aymmetric airfoils"]

for i, (s1_set, s2_set) in enumerate(zip([sym_s1, asym_s1], [sym_s2, asym_s2])):
    plt.figure(figsize=(7,3.5))
    for s1, s2 in zip(s1_set, s2_set):
        z_a = np.array([R(t,s1,s2) for t in theta])*np.exp(1j*theta) 
        l = np.real( zeta(z_a)[0] - zeta(z_a)[n_points//2])
        plt.plot(np.real(zeta(z_a))/l - np.real(zeta(z_a)[n_points//2])/l, np.imag(zeta(z_a))/l, label=f"s1/c = {round(s1,3)}, s2/c = {round(s2,3)}")

        h_x = np.array( [np.real(zeta(z_a)[i]) for i in range(n_points//2)])
        h_y = np.array( [np.imag(zeta(z_a)[i] + zeta(z_a)[n_points-1-i])/2 for i in range(n_points//2)])

    plt.plot(h_x/l-np.real(zeta(z_a)[n_points//2])/l, h_y/l, '--', color='black')
    plt.axis('equal')
    plt.xlabel('x/l')
    plt.ylabel('y/l')
    plt.title(Titles[i])
    plt.grid(True)
    plt.legend()
    plt.show()



def print_camber(transformation, z_arr):
    h_x = np.array( [np.real(transformation(z_arr)[i]) for i in range(n_points//2)])
    h_y = np.array( [np.imag(transformation(z_arr)[i] + transformation(z_arr)[n_points-1-i])/2 for i in range(n_points//2)])
    print(h_x/l-np.real(transformation(z_arr)[n_points//2])/l)
    print(h_y/l)