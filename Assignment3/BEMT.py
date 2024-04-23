import numpy as np

#Initial conditions

T = 10
rho = 1.225
A = 0.785
R = 0.5
N_b = 2
Omega = 100
theta = 5*np.pi/180
c = 0.1

C_d = 0.5
c_l = 0.01
v_i = np.sqrt(T/(2*rho*A))

y = R/2

for i in range(10):

    U_T = Omega * y
    U_P = v_i

    phi = U_P / U_T
    alpha = theta - phi

    r = y / R

    sigma = N_b * c /(np.pi * R)

    f = N_b / 2 * ((1 - r)/(r * phi))
    F = (2 / np.pi) * np.cos(np.exp(-f))**(-1)

    lamba = sigma * c_l / (16 * F) * (np.sqrt(1 + 32/(sigma * c_l) * theta * r) - 1)


    # C_T = 4* np.integrate(lamba **2 * r, r=0, r=1)
    # C_T = 4* np.integrate(lamba **3 * r, r=0, r=1)

    v_i = U_T * lamba

    print(lamba)