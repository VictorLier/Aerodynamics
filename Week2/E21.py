import sympy as sp

R, Pc, Pinf, rho, U, theta = sp.symbols('R Pc Pinf rho U theta')

pp = 1/2 * rho * U**2 * (1-4*sp.sin(theta)**2)

L = -sp.integrate(pp*sp.sin(theta)*R, (theta, 0, sp.pi))

print(L)