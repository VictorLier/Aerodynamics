import numpy as np
import math
import os
import subprocess

def solve_fourier_coefs(nA: int, m0: float, b: float, c: np.ndarray, alpha_L0: np.ndarray, alpha: np.ndarray, theta_array: np.ndarray) -> np.ndarray:
    '''
    Calculates the Fourier coefficients for a given airfoil at a specific angle of attack
    nA: Number of Fourier coefficients
    m0: Slope of lift curve
    b: Wing span
    c: 1D array of cord values
    alpha_L0: 1D array of zero lift angle of attack
    alpha: 1D array of angle of attack
    theta_array: 1D array of theta values
    '''
    # Initialize ordinate vector - o
    o = np.full(nA, (-alpha + alpha_L0) )

    # Initialize coefficient matrix - B
    B = np.zeros((nA, nA))

    # Fill in B matrix -  Silde 26 week 8
    for j, theta in enumerate(theta_array):
        for i in range(nA):
            B[j, i] = (-4*b*np.sin((1+i)*theta) / (m0 * c[j]) - (1+i)*(np.sin((1+i)*theta)/np.sin(theta))) # Silde 26 week 8

    # Solve for A
    A = np.linalg.solve(B, o)

    A[abs(A) < 1e-20] = 0
    return A


def get_xfoil_data(airfoil_name: str, AOA_start: float, AOA_end: float, AOA_step: float, Re: float):
    '''
    Run xfoil and get data for a given airfoil
    Input: airfoil_name, AOA_start [degrees], AOA_end [degrees], AOA_step [degrees], Re
    Output: alpha [radians], alpha0 [radians], m0 (slope), cd_friction
    '''
    # File name
    polar_file_path = f"Assignment2/{airfoil_name}_{AOA_start}_{AOA_end}_{AOA_step}_{Re}.txt"



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
    input_file = open("Assignment2/input.in", 'w')
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



class Battery:
    def __init__(self, mass = None, energy = None, power = None):
        '''
        Creates battery object with mass, energy, or power attributes.
        Rounds up the mass to meet the requirements of the power and energy attributes.
        mass: int, mass of battery in kg
        energy: minimum energy of battery in J
        power: minimum power of battery in W
        '''
        if sum(arg is not None for arg in [mass, energy, power]) != 1:  # Checks if to many attributes are given
            raise ValueError("Exactly one of mass, energy, or power must be specified")
        
        if mass is not None: # If mass is given, calculate energy and power
            if mass <= 0 or mass % 1 != 0: # Checks if mass is positive and an integer
                raise ValueError("Mass must be positive and an integer")
            self.mass = mass # kg
            self.energy = 750000 * self.mass # J
            self.power = 600 * self.mass # W

        elif energy is not None: # If energy is given, calculate mass and power
            if energy <= 0: # Checks if energy is positive
                raise ValueError("Energy must be positive")
            mass = energy / 750000
            self.mass = math.ceil(mass)
            self.energy = 750000 * self.mass # J
            self.power = 600 * self.mass # W

        elif power is not None: # If power is given, calculate mass and energy
            if power <= 0: # Checks if power is positive
                raise ValueError("Power must be positive")
            mass = power / 600
            self.mass = math.ceil(mass) # kg
            self.energy = 750000 * self.mass # J
            self.power = 600 * self.mass # W


class Prop:
    def __init__(self, radius: float, inner_radius: float, no_blades: int, tip_speed: float, cord: float, alpha: float, airfoil: str)-> None:
        '''
        radius: float, radius of the propeller in m
        inner_radius: float, inner radius of the propeller in m
        no_blades: int, number of blades on the propeller
        tip_speed: float, tip speed of the propeller in m/s
        cord: float, cord of the propeller in m
        alpha: float, angle of attack of the propeller in degrees
        airfoil: str, airfoil used on the propeller
        '''
        self.radius = radius
        self.inner_radius = inner_radius
        self.no_blades = no_blades
        self.tip_speed = tip_speed
        self.cord = cord
        self.alpha = alpha
        self.airfoil = airfoil
        self.rho = 1.225 # kg/m^3

    def ideal_thrust(self):
        '''
        Returns the ideal thrust of the propeller - Week 9 slide 8
        '''
        T = 2 * self.rho * self.area * self.tip_speed**2
        return T
    

    def ideal_power(self):
        '''
        Returns the ideal power of the propeller - Week 9 slide 8
        '''
        P = 2 * self.rho * self.area * self.tip_speed**3
        return P

    def get_clalpha(self):
        '''
        Returns the lift coefficient of the propeller
        '''
        self.clalpha = 2 * math.pi # Lift coefficient approx  --  Ved ikke med Reynolds tal
    
    def solidity(self):
        '''
        Returns the solidity of the propeller
        '''
        self.sigma = self.no_blades * self.cord / ( math.pi * self.radius)

    def BEMT(self, twist=0, n_elements=10, tip_loss=True):
        '''
        Calculates the thrust and power of the propeller using the Blade Element Momentum Theory
        twist: int, 0 = constant, 1 = linear, 2 = ideal
        n_elements: int, number of elements to divide the propeller into
        tip_loss: bool, if True, tip loss is accounted for
        '''
        self.solidity()
        self.get_clalpha()

        # Create arrays for the elements
        r = np.linspace(self.inner_radius, self.radius, n_elements)

        # Calculate the twist
        if twist == 0:
            theta = 1
        elif twist == 1:
            theta_0 = 0 # Root twist
            theta_tw = 2
            theta = theta_0 + r * theta_tw
        elif twist == 2:
            theta_tip = self.alpha
            theta = theta_tip / r
        else:
            raise ValueError("Invalid twist value")
    
        # Calculate lambda
        if tip_loss:
            F_init = 1
            F = 0
            while np.abs(np.max(np.abs(F_init) - np.abs(F))) > 0.0001:
                F = F_init
                lambda_r = self.sigma * self.clalpha / (16 * F) * (np.sqrt(1 + (32 * F) / (self.sigma * self.clalpha) * theta * r) - 1) # 3.126
                phi = lambda_r / r
                f = self.no_blades / 2 * ((1 - r) / (r * phi))
                F_init = (2 / np.pi) * np.cos(np.exp(-f))**(-1)
        else:
            lambda_r = self.sigma * self.clalpha / (16) * (np.sqrt(1 + (32) / (self.sigma * self.clalpha) * theta * r) - 1) # 3.126

        # Calculate the thrust and power
        self.CT = 4 * np.sum(lambda_r**2 * r) # Lidt et skud i tågen
        self.CP = 4 * np.sum(lambda_r**3 * r) # Lidt et skud i tågen
        

class Wing:
    def __init__(self, min_lift: float, min_speed: float, max_speed: float, center_chord: float, airfoil = "NACA4415"):
        '''
        min_lift: float, minimum lift of the wing in N
        min_speed: float, minimum speed of the wing in m/s
        max_speed: float, maximum speed of the wing in m/s
        chord: float, chord of the wing in m
        airfoil: str, airfoil of the wing eg. "NACA4415"
        '''
        self.min_lift = min_lift
        self.min_speed = min_speed
        self.max_speed = max_speed
        self.c0 = center_chord
        self.airfoil = airfoil
        self.rho = 1.225 # Air density [kg/m^3]

    def Re(self, speed: float):
        '''
        calculates the Reynolds number of the wing
        '''
        v = 1.5e-5 # Dynamic viscosity of air [m^2/s]
        self.Re = speed * self.c0 / v # Reynolds number
    
    def optimal_ClCd(self):
        '''
        Finds the optimal Cl/Cd ratio of the wing at min speed and the Cl and Cd values at this angle of attack
        '''
        self.Re()

    def dimensions(self):
        '''
        Calculates the dimensions of the wing
        '''
        self.ClCd()
        self.area = 2 * self.min_lift / (self.min_speed**2 * self.rho * self.Cl) # Wing area
        self.b = 4 * self.area / (np.pi * self.c0) # Wing span
        self.AR = 4 * self.b / (np.pi * self.c0) # Aspect ratio

    def max_speed_ClCd(self):
        '''
        Calculates the Cl and Cd values at max speed
        '''
        

    def analytical_CDi(self, speed=0):
        '''
        Calculates the induced drag of the wing
        '''

        self.ClCd()
        CDi = 1 / (np.pi * self.AR) * self.Cl**2



class Plane:
    def __init__(self, mass: float, battery: Battery, hover_prop: Prop, tilt_prop: Prop):
        self.mass = mass
        self.battery = battery
        self.hover_prop = hover_prop
        self.tilt_prop = tilt_prop
    
    def payload(self):
        payload = self.mass - self.battery.mass
        return payload

    def ideal_hover_time(self):
        '''
        Returns the ideal hover time with the power from - week 9 slide 8
        '''
        P_hover = self.battery.energy / self.hover_prop.ideal_power()
        return P_hover

    def required_hover_thrust(self):
        '''
        Returns the required thrust for the plane to hover
        '''
        T = self.mass * 9.81
        return T

if __name__ == "__main__":
    if True: # Test
        P1 = Prop(0.01, 0.2, 3, 100, 0.05, 1, "NACA 0012")
        P1.BEMT(2, 20)
        print(P1.CT)
        print(P1.CP)

        print("stop")
