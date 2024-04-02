import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

def solve_fourier_coefs(nA: int, m0: float, b: float, c: np.ndarray, alpha_L0: np.ndarray, alpha: np.ndarray, theta_array: np.ndarray) -> np.ndarray:
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
    # Initialize ordinate vector - o

    # b = 10

    o = np.full(nA, (alpha_L0 - alpha))

    # Initialize coefficient matrix - B
    B = np.zeros((nA, nA))

    # Fill in B matrix
    for j, theta in enumerate(theta_array):
        for i in range(nA):
            B[i, j] = (-4*b/(m0*c[j]) * np.sin(((1+i))*theta) - (i+1)*(np.sin(((1+i))*theta)/np.sin(theta))) # Silde 26 week 8

    # Solve for A
    A = np.linalg.solve(B, o)

    return A


def get_xfoil_data(airfoil_name: str, AOA_start: float, AOA_end: float, AOA_step: float, Re: float):
    '''
    Run xfoil and get data for a given airfoil
    Input: airfoil_name, AOA_start [degrees], AOA_end [degrees], AOA_step [degrees], Re
    Output: alpha [radians], alpha0 [radians], m0 (slope), cd_friction
    '''
    # File name
    polar_file_path = f"Assignment2\{airfoil_name}_{AOA_start}_{AOA_end}_{AOA_step}_{Re}.txt"



    # Check if polar file already exists
    if os.path.exists(polar_file_path):
        data = np.loadtxt(polar_file_path, skiprows=12)
        alpha = data[:, 0] * np.pi/180
        cl = data[:, 1]
        cd = data[:, 2]
        cdp = data[:, 3]

        # calculate slope and intercept of lift curve
        m0, cl0 = np.polyfit(alpha, cl, 1) # Assumes linear lift curve slope
        alpha0 = -cl0/m0 # Zero lift angle of attack in degrees
        cd_friction = cd # Friction drag coefficient

        return alpha, alpha0, m0, cd_friction
    
    # Create input file
    input_file = open("Assignment2\input.in", 'w')
    input_file.write(f"{airfoil_name}\n")
    input_file.write("PANE\n")
    input_file.write("OPER\n")
    input_file.write(f"Visc {Re}\n")
    input_file.write("PACC\n")
    input_file.write(f"{polar_file_path}\n\n")
    input_file.write("ITER 100\n")
    input_file.write("ASeq {0} {1} {2}\n".format(AOA_start, AOA_end, AOA_step))
    input_file.write("\n\n")
    input_file.write("quit\n")
    input_file.close()

    # Run xfoil
    subprocess.call("Assignment2\\xfoil.exe < Assignment2\\input.in", shell=True)

    # Read data from polar file
    data = np.loadtxt(polar_file_path, skiprows=12)
    alpha = data[:, 0] * np.pi/180
    cl = data[:, 1]
    cd = data[:, 2]
    cdp = data[:, 3]

    # calculate slope and intercept of lift curve
    m0, cl0 = np.polyfit(alpha, cl, 1) # Assumes linear lift curve slope
    alpha0 = -cl0/m0 # Zero lift angle of attack in degrees
    cd_friction = cd # Friction drag coefficient

    return alpha, alpha0, m0, cd_friction

class EllipticAirfoil():
    def __init__(self, AR, airfoil_name = "NACA4415", Re = 6e6, nA = 10, AOA_start = -6, AOA_end = 10, AOA_step = 0.5):
        '''
        airfoil_name: Airfoil section name eg. NACA4415
        Re: Reynolds number
        nA: Number of Fourier coefficients/theta values
        AR: Aspect ratio
        AOA_start: Start angle of attack [degrees]
        AOA_end: End angle of attack [degrees]
        AOA_step: Angle of attack step [degrees]
        '''
        self.airfoil_name = airfoil_name
        self.Re = Re
        self.nA = nA
        self.AR = AR
        self.AOA_start = AOA_start
        self.AOA_end = AOA_end
        self.AOA_step = AOA_step

        # Get xfoil data - list of angles [radians], zero lift angle of attack [radians], slope of lift curve, friction drag coefficient
        self.AOA, self.alpha0, self.m0, self.cd_friction = get_xfoil_data(self.airfoil_name, self.AOA_start, self.AOA_end, self.AOA_step, self.Re)

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

class ConstantAirfoil(EllipticAirfoil):
    def __init__(self, AR, airfoil_name = "NACA4415", Re = 6e6, nA = 10, AOA_start = -6, AOA_end = 10, AOA_step = 0.5):
        self.c = 1
        super().__init__(AR, airfoil_name, Re, nA, AOA_start, AOA_end, AOA_step)
        self.b = AR

    def ai(self) -> np.ndarray:
        '''
        Calculate the induced angle of attacks at each theta
        return: Induced angle of attack, Angle of attack in degrees
        '''
        ai = np.zeros((len(self.AOA), self.nA))
        for i, AOA in enumerate(self.AOA):
            for j, theta in enumerate(self.theta):
                for n in range(self.nA):
                    ai_sum =+ (n+1)*self.A[n,i] * np.sin((n+1)*theta) / np.sin(theta) # slide 26 week 8
                ai[i, j] = ai_sum
        return ai, self.AOA*180/np.pi



if __name__ == "__main__":

    if False: # Task 1
        AR = [4, 6, 8, 10, 10000]
        # AR = [10]
        Foils = []
        for ar in AR:
            Foils.append(EllipticAirfoil(AR = ar))
        
        # General Cl
        plt.figure()
        for foil in Foils:
            plt.plot(foil.general_Cl()[1], foil.general_Cl()[0])
        plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
        plt.xlabel("Angle of attack (degrees)")
        plt.ylabel("Lift coefficient")
        plt.title("Lift coefficient vs. Angle of attack for different aspect ratios")
        plt.grid(True)

        # General Cdi
        plt.figure()
        for foil in Foils:
            plt.plot(foil.general_Cdi()[1], foil.general_Cdi()[0])
        plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
        plt.xlabel("Angle of attack (degrees)")
        plt.ylabel("Induced drag coefficient")
        plt.title("Induced drag coefficient vs. Angle of attack for different aspect ratios")
        plt.grid(True)

        # ai
        plt.figure()
        for foil in Foils:
            plt.plot(foil.ai()[1], foil.ai()[0])
        plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
        plt.xlabel("Angle of attack (degrees)")
        plt.ylabel("Induced angle of attack")
        plt.title("Induced angle of attack vs. Angle of attack for different aspect ratios")
        plt.grid(True)

        # AOA
        plt.figure()
        for foil in Foils:
            plt.plot(foil.AOA*180/np.pi, foil.cd_friction)
        plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
        plt.xlabel("Angle of attack (degrees)")
        plt.ylabel("Friction drag coefficient")
        plt.title("Friction drag coefficient vs. Angle of attack for different aspect ratios")
        plt.grid(True)



        plt.show()

    if True: # Task 2
        if True: # Part A
            AR = [4, 6, 8, 10, 10000]
            # AR = [10]
            Foils = []
            for ar in AR:
                Foils.append(ConstantAirfoil(AR = ar, AOA_start=0, AOA_end=10, AOA_step=0.5))

            # AOA 0
            plt.figure()
            for foil in Foils:
                plt.plot(foil.theta, abs(foil.ai()[0][0]*180/np.pi))
            plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
            plt.xlabel("Theta")
            plt.ylabel("Induced angle of attack")
            plt.title("Induced angle of attack vs. Theta for different aspect ratios and angle of attack = 0 degrees")
            plt.grid(True)

            # AOA 5
            plt.figure()
            for foil in Foils:
                plt.plot(foil.theta, foil.ai()[0][1]*180/np.pi)
            plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
            plt.xlabel("Theta")
            plt.ylabel("Induced angle of attack")
            plt.title("Induced angle of attack vs. Theta for different aspect ratios and angle of attack = 5 degrees")
            plt.grid(True)

            # AOA 10
            plt.figure()
            for foil in Foils:
                plt.plot(foil.theta, foil.ai()[0][2]*180/np.pi)
            plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
            plt.xlabel("Theta")
            plt.ylabel("Induced angle of attack")
            plt.title("Induced angle of attack vs. Theta for different aspect ratios and angle of attack = 10 degrees")
            plt.grid(True)
            
            plt.show()
            
        if False: # Part B
            AR = [4, 6, 8, 10, 10000]
            Foils = []

            for ar in AR:
                Foils.append(ConstantAirfoil(AR = ar, AOA_start=-4))

            # CL
            plt.figure()
            for foil in Foils:
                plt.plot(foil.general_Cl()[1], foil.general_Cl()[0])
            plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
            plt.xlabel("Angle of attack (degrees)")
            plt.ylabel("Lift coefficient")
            plt.title("Lift coefficient vs. Angle of attack for different aspect ratios")
            plt.grid(True)

            # Cdi
            plt.figure()
            for foil in Foils:
                plt.plot(foil.general_Cdi()[1], foil.general_Cdi()[0])
            plt.legend(["AR=4", "AR=6", "AR=8", "AR=10", "AR=inf"])
            plt.xlabel("Angle of attack (degrees)")
            plt.ylabel("Induced drag coefficient")
            plt.title("Induced drag coefficient vs. Angle of attack for different aspect ratios")
            plt.grid(True)

            plt.show()
    

    if False: # test
        alpha, alpha0, m0, cd_friction = get_xfoil_data("NACA4415", -6, 10, 0.5, 6e6)
        


    print("Stop")