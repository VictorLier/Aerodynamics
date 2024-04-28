import numpy as np
import math
import os
import subprocess
import matplotlib.pyplot as plt
from scipy.stats import linregress


def linear_reg_limits(x: np.ndarray, y: np.ndarray, start: int, min_step: int, buffer: float=1):
    '''
    Finds the linear regression of a set of data points and the limits for a valid data set
    x: 1D array of x values
    y: 1D array of y values
    start: int, start index
    min_steps: int, first step size
    outputs
    slope: float, slope of the linear regression
    intercept: float, intercept of the linear regression
    min_step: int, end index
    buffer: float, amount of min_steps to go backwards
    '''

    stop = len(x[start:])

    i = 1
    r_square = 1

    while r_square > 0.99:
        x_data = x[start:start+min_step+i]
        y_data = y[start:start+min_step+i]
        slope, intercept, r_value, p_value, std_err = linregress(x_data, y_data) 
        r_square = r_value**2
        i += 1
        if start+min_step+i >= stop:
            break

    back = int(min_step*buffer)
    x_data = x[start:start+min_step+i-back]
    y_data = y[start:start+min_step+i-back]
    slope, intercept, r_value, p_value, std_err = linregress(x_data, y_data)
    max_step = start+min_step+i-back

    return slope, intercept, max_step

def xfoil(re: float, airfoil: str, AOA_start: float=-10, AOA_end: float=20, AOA_step: float=0.1):
    '''
    Runs xfoil to get foil data
    re: float, Reynolds number
    airfoil: str, airfoil name
    AOA_start: float, start angle of attack [deg]
    AOA_end: float, end angle of attack [deg]
    AOA_step: float, step size of angle of attack [deg]
    outputs
    alpha: 1D array of angle of attack [rad]
    cl: 1D array of lift coefficient
    cd: 1D array of drag coefficient
    m0: float, lift curve slope
    alpha0: float, zero lift angle of attack
    stall_index: int, index of stall angle of attack
    '''
    # Creates file name
    polar_file_path = f"Assignment3/{airfoil}_{AOA_start}_{AOA_end}_{AOA_step}_{re}.txt"

    # If polar file already exists
    if os.path.exists(polar_file_path):
        data = np.loadtxt(polar_file_path, skiprows=12)

    else: # File does not exist, xfoil needs to run
        # Create input file
        input_file = open("Assignment3/input.in", 'w')
        input_file.write(f"{airfoil}\n")
        input_file.write("PANE\n")
        input_file.write("OPER\n")
        input_file.write(f"Visc {re}\n")
        input_file.write("PACC\n")
        input_file.write(f"{polar_file_path}\n\n")
        input_file.write("ITER 100\n")
        input_file.write("ASeq {0} {1} {2}\n".format(AOA_start, AOA_end, AOA_step))
        input_file.write("\n\n")
        input_file.write("quit\n")
        input_file.close()
    
        # Run xfoil
        subprocess.call("Assignment3\\xfoil.exe < Assignment3\\input.in", shell=True)

    # Read data from polar file
    data = np.loadtxt(polar_file_path, skiprows=12)
    alpha = data[:, 0] * np.pi/180 # Angle of attack
    cl = data[:, 1] # Lift coefficient
    cd = data[:, 2] # Friction drag coefficient
    cdp = data[:, 3] # Pressure drag coefficient
    cdf = cd - cdp # Friction drag coefficient

    # Save data and calculate slope and intercept of lift curve
    # Lift curve slope and zero lift angle of attack - Assumes linear lift curve slope and Find stall AOA
    start = len(alpha)//6 # Start point of linear regression
    m0, cl0, stall_index = linear_reg_limits(alpha, cl, start, 5)
    alpha0 = -cl0/m0 # Zero lift angle of attack in degrees
    
    return alpha, cl, cdf, m0, alpha0, stall_index

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
    def __init__(self, radius: float, inner_radius: float = 0.02, no_blades: int = 2, tip_speed: float = 100, cord: float = 0.02, alpha: float = 2, airfoil: str = "NACA23012")-> None:
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
    def __init__(self, min_lift: float, min_speed: float = 11, max_speed: float = 22, aspect_ratio: float = 10, airfoil = "NACA4415",optimum_speed:float = 18):
        '''
        min_lift: float, minimum lift of the wing in Kg
        min_speed: float, minimum speed of the wing in m/s
        max_speed: float, maximum speed of the wing in m/s
        aspect_ratio: float, aspect ratio of the wing
        airfoil: str, airfoil of the wing eg. "NACA4415"
        optimum_speed: float, optimal speed of the wing in m/s
        bodyarea: float, body area in m^2
        bodyCD: float, body drag coefficient
        '''
        self.min_lift = min_lift * 9.82 # Minimum lift in N
        self.min_speed = min_speed # Minimum speed in m/s
        self.max_speed = max_speed # Maximum speed in m/s
        self.AR = aspect_ratio # Aspect ratio
        self.airfoil = airfoil # Airfoil of the wing
        self.rho = 1.225 # Air density [kg/m^3]
        self.optimum_speed = optimum_speed # Optimal speed in m/s

    def reynolds(self, speed: float, c = 0.6): # 
        '''
        calculates the Reynolds number of the wing
        speed: float, speed of the wing in m/s
        c: float, chord of the wing in m - Assumes 0.6 as default
        '''
        v = 1.5e-5 # Dynamic viscosity of air [m^2/s]
        re = speed * c / v # Reynolds number
        return re

    def xfoil_data(self, plot=False):
        '''
        Runs xfoil and gets data for the wing
        plot: bool, if True, plots the data
        '''
        # Get Reynolds number for xfoil
        re = self.reynolds(self.optimum_speed)
        
        # Run xfoil
        self.alpha, self.cl, self.cd, self.m0, self.alpha0, self.stall_index = xfoil(re, self.airfoil)

        if plot:
            plt.figure()
            plt.title("Xfoil data")
            plt.plot(self.alpha * 180/np.pi, self.cl, label="Cl")
            plt.plot(self.alpha[self.stall_index] * 180/np.pi, self.cl[self.stall_index], 'ro', label="Stall angle of attack: {:.2f}".format(self.alpha[self.stall_index] * 180/np.pi))
            plt.xlabel("Angle of attack [deg]")
            plt.ylabel("Coefficient")
            plt.legend()
            plt.grid()

            plt.figure()
            plt.title("Xfoil data")
            plt.plot(self.alpha * 180/np.pi, self.cd, label="Cd")
            plt.xlabel("Angle of attack [deg]")
            plt.ylabel("Coefficient")
            plt.legend()
            plt.grid()
        
            plt.figure()
            plt.title("Xfoil data")
            plt.plot(self.alpha * 180/np.pi, self.cl / self.cd, label="Cl/Cd")
            plt.plot(self.alpha[self.stall_index] * 180/np.pi, self.cl[self.stall_index] / self.cd[self.stall_index], 'ro', label="Stall angle of attack: {:.2f}".format(self.alpha[self.stall_index] * 180/np.pi))
            plt.xlabel("Angle of attack [deg]")
            plt.ylabel("Coefficient")
            plt.legend()
            plt.grid()

    def coefficient_lift(self):
        '''
        Calculates the lift coefficient of the wing using lifting line theory
        '''
        # get data from xfoil
        self.xfoil_data()
        self.CL = np.pi * self.m0 * (self.alpha - self.alpha0) * self.AR / (self.AR * np.pi + self.m0) # Lift coefficient
    
    def coefficient_drag(self):
        '''
        Calculates the drag coefficient of the wing using lifting line theory
        '''
        self.coefficient_lift()
        CDi = 1 / (np.pi * self.AR) * self.CL**2 # Induced drag coefficient
        self.CD = CDi + self.cd # Total drag coefficient

    def optimal_ClCd(self, plot=False):
        '''
        Finds the optimal Cl/Cd ratio of the wing at min speed and the Cl and Cd values at this angle of attack
        '''
        self.coefficient_drag()

        CLCD = self.CL / self.CD

        self.index_optimal = np.argmax(CLCD)
        self.CL_optimal = self.CL[self.index_optimal]
        self.CD_optimal = self.CD[self.index_optimal]
        self.alpha_optimal = self.alpha[self.index_optimal]


        if plot:
            plt.figure()
            plt.title("Optimal Cl/Cd")
            plt.plot(self.alpha * 180/np.pi, self.CL, label="CL")
            plt.plot(self.alpha * 180/np.pi, self.CD, label="CD")
            plt.xlabel("Angle of attack [deg]")
            plt.ylabel("Coefficient")
            plt.legend()
            plt.grid()

            plt.figure()
            plt.title("Optimal Cl/Cd")
            plt.plot(self.alpha * 180/np.pi, CLCD, label="Cl/Cd")
            plt.plot(self.alpha_optimal * 180/np.pi, self.CL_optimal / self.CD_optimal, 'ro', label="Optimal Cl/Cd - angle of attack: {:.2f}".format(self.alpha_optimal * 180/np.pi))
            plt.xlabel("Angle of attack [deg]")
            plt.ylabel("Coefficient")
            plt.legend()
            plt.grid()

    def dimensions(self, plot = False):
        '''
        Calculates the dimensions of the wing such that the optimal CD/Cd ratio is achieved at the optimum speed
        '''
        self.optimal_ClCd()
        self.S = 2 * self.min_lift / (self.optimum_speed**2 * self.rho * self.CL_optimal) # Wing area
        self.b = np.sqrt(self.AR * self.S) # Wing span
        self.c0 = 4 * self.S / (np.pi * self.b) # Root chord

        # Check if minimum speed is achievable
        CL_required = self.min_lift / (0.5 * self.rho * self.min_speed**2 * self.S)
        if CL_required > self.CL[self.stall_index]:
            raise ValueError("Minimum speed not achievable")

        if plot: # plot wing shape
            angles = np.linspace(0, 2*np.pi, 100)
            x = self.b * np.cos(angles)
            y = self.c0 * np.sin(angles)
            plt.figure()
            plt.plot(x, y)
            plt.title("Wing shape")
            plt.xlabel("Span [m]")
            plt.ylabel("Chord [m]")
            plt.axis('equal')
            plt.grid()

    def maxspeed_ClCd(self, plot=False):
        '''
        Calculates the Cl, aoa an CD of the wing at max speed
        '''
        self.dimensions()
        self.CL_max = self.min_lift / (0.5 * self.rho * self.max_speed**2 * self.S) # Lift coefficient
        self.index_max = np.argmin(np.abs(self.CL - self.CL_max))
        self.CD_max = self.CD[self.index_max]
        self.alpha_max = self.alpha[self.index_max]

        if plot:
            plt.plot(self.alpha_max * 180/np.pi, self.CL[self.index_max] / self.CD[self.index_max], 'ro', label="Max speed Cl/Cd - angle of attack: {:.2f}".format(self.alpha_max * 180/np.pi))
            plt.legend()

    def minspeed_ClCd(self, plot=False):
        '''
        Calculates the Cl, aoa an CD of the wing at min speed
        '''
        self.dimensions()
        self.CL_min = self.min_lift / (0.5 * self.rho * self.min_speed**2 * self.S) # Lift coefficient
        self.index_min = np.argmin(np.abs(self.CL - self.CL_min))
        self.CD_min = self.CD[self.index_min]
        self.alpha_min = self.alpha[self.index_min]

        if plot:
            plt.plot(self.alpha_min * 180/np.pi, self.CL[self.index_min] / self.CD[self.index_min], 'ro', label="Min speed Cl/Cd - angle of attack: {:.2f}".format(self.alpha_min * 180/np.pi))
            plt.legend()

    def CLCD_speed(self, plot=False, speed: np.ndarray = None):
        '''
        Calculates the Cl, aoa an CD of the wing at different speeds
        speed: 1D array of speeds to calculate the Cl, aoa and CD at. Default is from max speed to 0.1 m/s
        '''
        self.dimensions()

        #Initialize arrays
        if speed is None:
            speed_ar = np.linspace(0.1, self.max_speed, 100)
            self.CL_speed = np.zeros(len(speed_ar))
            self.CD_speed = np.zeros(len(speed_ar))
            self.alpha_speed = np.zeros(len(speed_ar))
        else:
            speed_ar = speed
            self.CL_speed = np.zeros(len(speed_ar))
            self.CD_speed = np.zeros(len(speed_ar))
            self.alpha_speed = np.zeros(len(speed_ar))

        for i, speed in enumerate(speed_ar):
            CL_required = self.min_lift / (0.5 * self.rho * speed**2 * self.S)
            closest_index = np.argmin(np.abs(self.CL - CL_required))
            if closest_index < self.stall_index: # Check if stalling
                self.CL_speed[i] = self.CL[closest_index]
                self.CD_speed[i] = self.CD[closest_index]
                self.alpha_speed[i] = self.alpha[closest_index]
            else:
                self.CL_speed[i] = self.CL[self.stall_index]
                self.CD_speed[i] = self.CD[self.stall_index]
                self.alpha_speed[i] = self.alpha[self.stall_index]
        
        if plot:
            plt.figure()
            plt.title("Cl at different speeds")
            plt.plot(speed_ar, self.CL_speed, label="CL")
            plt.xlabel("Speed [m/s]")
            plt.ylabel("Coefficient")
            plt.legend()
            plt.grid()

            plt.figure()
            plt.title("Cd at different speeds")
            plt.plot(speed_ar, self.CD_speed, label="CD")
            plt.xlabel("Speed [m/s]")
            plt.ylabel("Coefficient")
            plt.legend()
            plt.grid()

            plt.figure()
            plt.title("Angle of attack at different speeds")
            plt.plot(speed_ar, self.alpha_speed * 180/np.pi, label="Angle of attack")
            plt.xlabel("Speed [m/s]")
            plt.ylabel("Angle of attack [deg]")
            plt.legend()
            plt.grid()

class Plane:
    def __init__(self, battery: Battery = Battery(4), center_rotor: Prop = Prop(0.3), tip_rotor: Prop = Prop(0.15), wing: Wing = None):
        self.battery = battery
        self.mass = self.battery.mass + 14 # kg
        self.center_rotor = center_rotor
        self.tip_rotor = tip_rotor
        if wing is None:
            self.wing = Wing(self.mass)
        else:
            self.wing = wing

        # Estimate of the body specifications
        self.body_area = self.center_rotor.radius * 2 + self.tip_rotor.radius * 0.5 * 0.1 # m
        self.body_CD = 0.2 # Drag coefficient

    def drag_force(self, plot=False):
        '''
        Calculates the drag force of the plane at different speeds
        '''
        # Create speed array
        speed = np.linspace(0.1, 22, 100)
        
        # Calculate wing CD at different speeds
        self.wing.CLCD_speed(speed=speed)

        # Calculate the drag force
        wing_drag = 0.5 * self.wing.rho * speed**2 * self.wing.S * self.wing.CD_speed

        # Calculate body drag at different speeds
        body_drag = 0.5 * self.wing.rho * speed**2 * self.body_area * self.body_CD

        # Total drag force
        total_drag = wing_drag + body_drag

        if plot:
            plt.figure()
            plt.title("Drag force")
            plt.plot(speed, wing_drag, label="Wing drag")
            plt.plot(speed, body_drag, label="Body drag")
            plt.plot(speed, total_drag, label="Total drag")
            plt.xlabel("Speed [m/s]")
            plt.ylabel("Force [N]")
            plt.legend()
            plt.grid()

    def lift_force(self, plot=False):
        '''
        Calculates the lift force of the plane at different speeds
        '''
        # Create speed array
        speed = np.linspace(0.1, 22, 100)

        # Calculate wing CL at different speeds
        self.wing.CLCD_speed(speed=speed)

        # Calculate the lift force
        wing_lift = 0.5 * self.wing.rho * speed**2 * self.wing.S * self.wing.CL_speed

        # Calculate the required center rotor thrust
        self.required_thrust = self.wing.min_lift - wing_lift

        if plot:
            plt.figure()
            plt.title("Lift force")
            plt.plot(speed, wing_lift, label="Wing lift")
            plt.plot(speed, self.required_thrust, label="Required thrust")
            plt.xlabel("Speed [m/s]")
            plt.ylabel("Force [N]")
            plt.legend()
            plt.grid()

    def fixed_aoa(self, alpha: float = -2, plot=False):
        '''
        Calculates the forces in x,y direction for a fixed angle of attack
        '''
        # Initilize wing parameters
        self.wing.dimensions()
        # Create speed array
        speed = np.linspace(0.1, 22, 100)

        # Finds index closest to the fixed angle of attack
        index = np.argmin(np.abs(self.wing.alpha - alpha*np.pi/180))

        # Calculate the lift force
        wing_lift = 0.5 * self.wing.rho * speed**2 * self.wing.S * self.wing.CL[index]

        # Calculate the drag force of the wing and body
        wing_drag = 0.5 * self.wing.rho * speed**2 * self.wing.S * self.wing.CD[index]
        body_drag = 0.5 * self.wing.rho * speed**2 * self.body_area * self.body_CD
        drag = wing_drag + body_drag

        # required thrust in x and y direction
        x_force = drag
        y_force = self.wing.min_lift - wing_lift

        # comined force
        force = np.sqrt(x_force**2 + y_force**2)

        if plot:
            plt.figure()
            plt.title("Forces at fixed angle of attack")
            plt.plot(speed, x_force, label="X-force")
            plt.plot(speed, y_force, label="Y-force")
            plt.plot(speed, force, label="Combined force")
            plt.xlabel("Speed [m/s]")
            plt.ylabel("Force [N]")
            plt.legend()
            plt.grid()


    def power(self):
        pass

if __name__ == "__main__":
    if True: # Test
        plane = Plane()
        plane.fixed_aoa(plot=True)

        plt.show()
    if False: # Question 3
        pass


    if False: # Question 4
        wing = Wing(min_lift=18, min_speed=11, max_speed=20, aspect_ratio=10)
        wing.xfoil_data(plot=True)
        wing.optimal_ClCd(plot=True)
        wing.maxspeed_ClCd(plot=True)
        wing.minspeed_ClCd(plot=True)
        print("optimal angle of attack: {:.2f}".format(wing.alpha_optimal * 180/np.pi))
        print("optimal Cl/Cd: {:.2f}".format(wing.CL_optimal / wing.CD_optimal))
        print("optimal Cl: {:.2f}".format(wing.CL_optimal))
        print("optimal Cd: {:.2f}".format(wing.CD_optimal))
        print("Max speed angle of attack: {:.2f}".format(wing.alpha_max * 180/np.pi))
        print("Max speed Cl: {:.2f}".format(wing.CL_max))
        print("Max speed Cd: {:.2f}".format(wing.CD_max))

        plt.show()

        print("stop")

    if False: # Question 5
        wing = Wing(min_lift=18, min_speed=11, max_speed=20, aspect_ratio=10)
        wing.CLCD_speed()

        plt.show()