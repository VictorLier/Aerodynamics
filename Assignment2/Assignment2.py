import numpy as np
import matplotlib.pyplot as plt

def solve_fourier_coefs(nA: int, m0: np.ndarray, b: float, c: np.ndarray, alpha_L0: np.ndarray, alpha: np.ndarray, theta_array: np.ndarray) -> np.ndarray:
    '''
    Calculates the Fourier coefficients for a given airfoil at a specific angle of attack
    nA: Number of Fourier coefficients
    m0: 1D array of slope of lift curve
    b: Wing span
    c: 1D array of cord values
    alpha_L0: 1D array of zero lift angle of attack
    alpha: 1D array of angle of attack
    theta: 1D array of theta values
    '''
    # Initialize ordinate vector - b
    o = (alpha_L0 - alpha)

    # Initialize coefficient matrix - B
    B = np.zeros((nA, nA))

    # Fill in B matrix
    for j, theta in enumerate(theta_array):
        for i in range(nA):
            B[i, j] = (-4*b/(m0[j]*c[j]) * np.sin(((1+i))*theta) - (np.sin(((1+i))*theta)/np.sin(theta)))

    # Solve for A
    A = np.linalg.solve(B, o)

    return A

class EllipticAirfoil():
    def __init__(self, airfoil_name: str, Re: float, nA: int, AR: int, AOA_int: int):
        '''
        airfoil_name: Airfoil section name eg. NACA4415
        Re: Reynolds number
        nA: Number of Fourier coefficients/theta values
        AR: Aspect ratio
        AOA_int: number of angle of attacks
        '''
        self.airfoil_name = airfoil_name
        self.Re = Re
        self.nA = nA
        self.AR = AR
        self.AOA_int = AOA_int

        # Zero lift angle of attack from xfoil - Should probably be automated
        self.alpha0 =  np.full(nA, -4.25*np.pi/180)

        # Slope of lift curve from xfoil - Should probably be automated
        self.m0 = np.full(nA, 2*np.pi)

        # List of angles of attack
        self.AOA = np.linspace(-6*np.pi/180, 10*np.pi/180, self.AOA_int)

        # List of theta values
        dpi = np.pi/self.nA
        self.theta = np.linspace(dpi, np.pi-dpi, self.nA)

        # Cord function - root cord c0 = 1
        c0 = 1
        self.c = np.sin(self.theta) * c0

        # Wing span
        self.b = self.AR * np.pi * c0 / 4

        # Wing area
        self.S = np.pi * c0/2 * self.b/2
        
        # Calculate Fourier coefficients
        # initialize A for each AOA
        self.A = np.empty((self.nA, len(self.AOA)))

        # Calculate A vector for each AOA and store in A
        for i, AOA in enumerate(self.AOA):
            self.A[:,i] = solve_fourier_coefs(nA=self.nA, m0=self.m0, b=self.b, c=self.c, alpha_L0=self.alpha0, alpha=AOA, theta_array=self.theta)

    def general_Cl(self) -> np.ndarray:
        '''
        Calculate the lift coefficients
        return: Lift coefficients, Angle of attacks in degrees
        '''
        # Calculate lift coefficient
        cl = np.pi * self.AR * self.A[0,:] # Slide 27 week 8
        return cl, self.AOA*180/np.pi
    
    def general_Cdi(self) -> np.ndarray:
        '''
        Calculate the induced drag coefficients
        return: Induced drag coefficients, Angle of attacks in degrees
        '''
        # Calculate lift coefficient - Slide 27 week 8
        for i in range(len(self.A)):
            cDi_sum =+ (i+1) * self.A[i,:]**2
        cDi = np.pi * self.AR * cDi_sum
        return cDi, self.AOA*180/np.pi

    def ai(self) -> np.ndarray:
        '''
        Calculate the maximum induced angle of attacks
        return: Induced angles of attack, Angle of attack in degrees
        '''
        for i in range(len(self.A)):    # Slide 26 week 8 - ??????????
            ai_sum =+ (i+1) * self.A[i,:]
        return ai_sum, self.AOA*180/np.pi


if __name__ == "__main__":

    eliptic4 = EllipticAirfoil("NACA4415", Re=1e6, nA=10, AR=4, AOA_int=100)
    eliptic6 = EllipticAirfoil("NACA4415", Re=1e6, nA=10, AR=6, AOA_int=100)
    eliptic8 = EllipticAirfoil("NACA4415", Re=1e6, nA=10, AR=8, AOA_int=100)
    eliptic10 = EllipticAirfoil("NACA4415", Re=1e6, nA=10, AR=10, AOA_int=100)
    elipticinf = EllipticAirfoil("NACA4415", Re=1e6, nA=10, AR=10000, AOA_int=100)

    plt.figure()
    plt.plot(eliptic4.general_Cl()[1], eliptic4.general_Cl()[0])
    plt.plot(eliptic6.general_Cl()[1], eliptic6.general_Cl()[0])
    plt.plot(eliptic8.general_Cl()[1], eliptic8.general_Cl()[0])
    plt.plot(eliptic10.general_Cl()[1], eliptic10.general_Cl()[0])
    plt.plot(elipticinf.general_Cl()[1], elipticinf.general_Cl()[0])
    
    plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
    plt.xlabel("Angle of attack (degrees)")
    plt.ylabel("Lift coefficient")
    plt.title("Lift coefficient vs. Angle of attack for different aspect ratios")
    plt.grid(True)
    plt.show()

    plt.figure()
    plt.plot(eliptic4.general_Cdi()[1], eliptic4.general_Cdi()[0])
    plt.plot(eliptic6.general_Cdi()[1], eliptic6.general_Cdi()[0])
    plt.plot(eliptic8.general_Cdi()[1], eliptic8.general_Cdi()[0])
    plt.plot(eliptic10.general_Cdi()[1], eliptic10.general_Cdi()[0])
    plt.plot(elipticinf.general_Cdi()[1], elipticinf.general_Cdi()[0])
    plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
    plt.xlabel("Angle of attack (degrees)")
    plt.ylabel("Induced drag coefficient")
    plt.title("Induced drag coefficient vs. Angle of attack for different aspect ratios")
    plt.grid(True)
    plt.show()

    plt.figure()
    plt.plot(eliptic4.ai()[1], eliptic4.ai()[0])
    plt.plot(eliptic6.ai()[1], eliptic6.ai()[0])
    plt.plot(eliptic8.ai()[1], eliptic8.ai()[0])
    plt.plot(eliptic10.ai()[1], eliptic10.ai()[0])
    plt.plot(elipticinf.ai()[1], elipticinf.ai()[0])
    plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
    plt.xlabel("Angle of attack (degrees)")
    plt.ylabel("Induced angle of attack")
    plt.title("Induced angle of attack vs. Angle of attack for different aspect ratios")
    plt.grid(True)
    plt.show()
    

    print("Stop")