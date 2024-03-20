def solve_for_A(theta_arr, m0=1, b=1, c=1, alpha_L0=1, alpha=1/2):
    B_temp = []
    for theta in theta_arr:
        B_temp.append([ (-4*b/(m0*c) * np.sin((1+i*2)*theta) - (np.sin((1+i*2)*theta)/np.sin(theta))) for i in range(len(theta_arr))])
    B = np.array(B_temp)
    b = np.ones(len(theta_arr))*(alpha_L0-alpha)
    a = np.linalg.inv(B) @ b
    return a

theta = np.linspace(np.pi/6, np.pi*5/6, 3)
print(solve_for_A(theta))