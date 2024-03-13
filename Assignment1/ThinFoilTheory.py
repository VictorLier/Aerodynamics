import matplotlib.pyplot as plt
import numpy as np

def Joukowski_CL(AoA, s1, s2, n=1.9):
    '''
    This function calculates the lift coefficient for a Joukowski airfoil
    '''
    cl = 4 * np.pi / n * np.sqrt(1 - 2 * s1 + s2**2) * (1 - ((1 + s1) / s1)**(-n)) * np.sin(AoA + np.arctan(s2 / (1 - s1)))
    return cl

def cl_symetric(AoA):
    '''
    Cl for a symmetric airfoil
    '''
    cl = 2 * np.pi * AoA
    return cl


def xy_to_xtheta(airfoil):
    '''
    This function takes in an airfoil in the form of x and y coordinates
    and returns the airfoil in the form of x, y and theta
    '''
    theta = np.arctan2(airfoil[1], airfoil[0] - 0.5)
    theta = np.flip(theta)
    airfoil = np.array([airfoil[0], airfoil[1], theta])
    return airfoil


def differtiate(airfoil):
    '''
    This function takes in an airfoil in the form of x, y and theta and returns 
    the airfoil in the form of x, y, theta and dy/dx
    '''
    dy_dx = np.gradient(airfoil[1], airfoil[0])
    airfoil = np.array([airfoil[0], airfoil[1], airfoil[2], dy_dx])
    return airfoil


def cl_asymetric(AoA,airfoil):
    '''
    this function takes in an airfoil in the form of x, y, theta and dy/dx and calculates the cl
    '''

    int = np.trapz(airfoil[3], x = airfoil[2])
    #int = np.trapz(airfoil[3])
    int2 = np.trapz(airfoil[3] * np.cos(airfoil[2]), x = airfoil[2])
    #int2 = np.trapz(airfoil[3] * np.cos(airfoil[2]))

    A0 = AoA - 1 / np.pi * int
    A1 = 2 / np.pi * int2

    cl = 2 * np.pi * (A0 + A1/2)
    return cl


if __name__ == "__main__":

    E_foil = np.array([[ 1.00000000e+00,  9.99083763e-01,  9.96346561e-01,  9.91809670e-01,
        9.85499105e-01,  9.77444881e-01,  9.67680337e-01,  9.56241541e-01,
        9.43166818e-01,  9.28496417e-01,  9.12272330e-01,  8.94538284e-01,
        8.75339908e-01,  8.54725056e-01,  8.32744293e-01,  8.09451511e-01,
        7.84904650e-01,  7.59166482e-01,  7.32305434e-01,  7.04396395e-01,
        6.75521465e-01,  6.45770612e-01,  6.15242191e-01,  5.84043289e-01,
        5.52289877e-01,  5.20106738e-01,  4.87627156e-01,  4.54992377e-01,
        4.22350819e-01,  3.89857064e-01,  3.57670636e-01,  3.25954600e-01,
        2.94873994e-01,  2.64594151e-01,  2.35278927e-01,  2.07088887e-01,
        1.80179491e-01,  1.54699308e-01,  1.30788311e-01,  1.08576286e-01,
        8.81813840e-02,  6.97088626e-02,  5.32500246e-02,  3.88813920e-02,
        2.66641228e-02,  1.66436851e-02,  8.84979053e-03,  3.29658612e-03,
       -1.69042110e-05, -1.10610470e-03],
       [9.85856858e-33, 3.29024414e-04, 1.31122986e-03, 2.93207934e-03,
       5.16756508e-03, 7.98453682e-03, 1.13411564e-02, 1.51874738e-02,
       1.94661179e-02, 2.41130942e-02, 2.90586821e-02, 3.42284180e-02,
       3.95441568e-02, 4.49251963e-02, 5.02894507e-02, 5.55546608e-02,
       6.06396215e-02, 6.54654121e-02, 6.99566119e-02, 7.40424812e-02,
       7.76580912e-02, 8.07453842e-02, 8.32541450e-02, 8.51428669e-02,
       8.63794961e-02, 8.69420386e-02, 8.68190173e-02, 8.60097664e-02,
       8.45245568e-02, 8.23845437e-02, 7.96215343e-02, 7.62775760e-02,
       7.24043678e-02, 6.80625004e-02, 6.33205368e-02, 5.82539453e-02,
       5.29439014e-02, 4.74759782e-02, 4.19387460e-02, 3.64223072e-02,
       3.10167891e-02, 2.58108248e-02, 2.08900469e-02, 1.63356236e-02,
       1.22228636e-02, 8.61991504e-03, 5.58658432e-03, 3.17329613e-03,
       1.42021416e-03, 3.56539529e-04]])

    E_foil = np.flip(E_foil, axis = 1)

    E_foil = xy_to_xtheta(E_foil)

    E_foil = differtiate(E_foil)

    AoA = np.linspace(-np.pi/18, np.pi/18, 100)
    AoA_degrees = np.degrees(AoA)

    E_foil = cl_asymetric(AoA, E_foil)

    Sym_cl = cl_symetric(AoA)

    JouKt_Analytical_cl = Joukowski_CL(AoA, 0.181, 0.260)

    plt.figure()
    plt.plot(AoA_degrees, Sym_cl, label = "Thin-airfoil symmetric")
    plt.plot(AoA_degrees, JouKt_Analytical_cl, label = "Analytical Cl")
    plt.plot(AoA_degrees, E_foil, label = "Thin-airfoil cambered")
    plt.xlabel('Angle of Attack [Degrees]')
    plt.ylabel('Coefficient of Lift (Cl) [-]')
    plt.legend()
    plt.grid()
    plt.show()

    print("stop")
