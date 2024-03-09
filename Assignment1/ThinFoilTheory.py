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

    test = np.array([[ 2.40000000e+00,  2.39541780e+00,  2.38172874e+00,  2.35903924e+00,
        2.32747940e+00,  2.28719932e+00,  2.23836574e+00,  2.18115904e+00,
        2.11577089e+00,  2.04240257e+00,  1.96126409e+00,  1.87257414e+00,
        1.77656091e+00,  1.67346372e+00,  1.56353546e+00,  1.44704565e+00,
        1.32428404e+00,  1.19556458e+00,  1.06122947e+00,  9.21653237e-01,
        7.77246475e-01,  6.28459124e-01,  4.75783066e-01,  3.19753859e-01,
        1.60951488e-01, -4.14130638e-16, -1.62434027e-01, -3.25644215e-01,
       -4.88888307e-01, -6.51393221e-01, -8.12361153e-01, -9.70976605e-01,
       -1.12641420e+00, -1.27784709e+00, -1.42445581e+00, -1.56543736e+00,
       -1.70001427e+00, -1.82744352e+00, -1.94702509e+00, -2.05810992e+00,
       -2.16010711e+00, -2.25249026e+00, -2.33480276e+00, -2.40666190e+00,
       -2.46776183e+00, -2.51787516e+00, -2.55685330e+00, -2.58462550e+00,
       -2.60119664e+00, -2.60664385e+00],
       [4.93038066e-32, 1.64548798e-03, 6.55760753e-03, 1.46636575e-02,
       2.58435722e-02, 3.99315637e-02, 5.67183947e-02, 7.59542592e-02,
       9.73522375e-02, 1.20592287e-01, 1.45325726e-01, 1.71180155e-01,
       1.97764761e-01, 2.24675942e-01, 2.51503180e-01, 2.77835086e-01,
       3.03265545e-01, 3.27399865e-01, 3.49860858e-01, 3.70294748e-01,
       3.88376819e-01, 4.03816718e-01, 4.16363312e-01, 4.25809022e-01,
       4.31993543e-01, 4.34806881e-01, 4.34191638e-01, 4.30144483e-01,
       4.22716784e-01, 4.12014338e-01, 3.98196218e-01, 3.81472708e-01,
       3.62102360e-01, 3.40388194e-01, 3.16673103e-01, 2.91334511e-01,
       2.64778386e-01, 2.37432689e-01, 2.09740370e-01, 1.82152041e-01,
       1.55118439e-01, 1.29082828e-01, 1.04473466e-01, 8.16962850e-02,
       6.11279112e-02, 4.31091614e-02, 2.79391344e-02, 1.58700097e-02,
       7.10265020e-03, 1.78309415e-03]])


    KT = np.array([[1.06355977, 1.05282914, 1.02120461, 0.96977022, 0.89989709, 0.81340948, 0.7129068, 0.60207506, 0.48580808, 0.37002603, 0.26118807, 0.16559884, 0.08867254, 0.03433098, 0.00467578],
                   [5.33394373e-33, 4.94451814e-03, 1.89531036e-02, 3.96821823e-02, 6.36456423e-02, 8.67826863e-02, 1.05133073e-01, 1.15512613e-01, 1.16074319e-01, 1.06653294e-01, 8.88262505e-02, 6.56644037e-02, 4.12130117e-02, 1.97826446e-02, 5.17694204e-03]]) 
    Jou = np.array([[1., 0.98804461, 0.95584563, 0.90559591, 0.83878413, 0.75712432, 0.66307107, 0.56011182, 0.45278599, 0.34642267, 0.24665724, 0.15884464, 0.08751551, 0.03601793, 0.00645877],
                    [-4.51529242e-32, 5.14471086e-03, 1.82365466e-02, 3.63973477e-02, 5.65011265e-02, 7.53705080e-02, 9.01024504e-02, 9.84020914e-02, 9.88937496e-02, 9.13607913e-02, 7.68545124e-02, 5.76201415e-02, 3.68168773e-02, 1.80555500e-02, 4.83620533e-03]])

    KT = np.flip(KT, axis = 1)
    Jou = np.flip(Jou, axis = 1)
    test = np.flip(test, axis = 1)

    KT = xy_to_xtheta(KT)
    Jou = xy_to_xtheta(Jou)
    test = xy_to_xtheta(test)

    KT = differtiate(KT)
    Jou = differtiate(Jou)
    test = differtiate(test)

    AoA = np.linspace(-np.pi/18, np.pi/18, 100)
    AoA_degrees = np.degrees(AoA)

    KT_cl = cl_asymetric(AoA, KT)
    Jou_cl = cl_asymetric(AoA, Jou)
    test_cl = cl_asymetric(AoA, test)

    Sym_cl = cl_symetric(AoA)

    JouKt_Analytical_cl = Joukowski_CL(AoA, 0.181, 0.260)

    plt.figure()
    #plt.plot(AoA_degrees, KT_cl, label = "KT")
    #plt.plot(AoA_degrees, Jou_cl, label = "Jou")
    plt.plot(AoA_degrees, Sym_cl, label = "Thin-airfoil symmetric")
    plt.plot(AoA_degrees, JouKt_Analytical_cl, label = "Analytical Cl")
    plt.plot(AoA_degrees, test_cl, label = "Thin-airfoil cambered")
    plt.legend()
    plt.show()

    print("stop")
