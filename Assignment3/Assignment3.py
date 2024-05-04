import numpy as np
import math
import os
import subprocess
import matplotlib.pyplot as plt
from scipy.stats import linregress

def linear_reg_limits(x: np.ndarray, y: np.ndarray, start: int, min_step: int, buffer: float=5):
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

    while r_square > 0.95:
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

def xfoil(re: float, airfoil: str, AOA_start: float=-10, AOA_end: float=20, AOA_step: float=0.1, plot=False):
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
    start = len(alpha)//8 # Start point of linear regression
    m0, cl0, stall_index = linear_reg_limits(alpha, cl, start, 10)
    alpha0 = -cl0/m0 # Zero lift angle of attack in degrees
    
    if plot:
        plt.figure()
        plt.plot(alpha * 180/np.pi, cl, label="Cl")
        plt.plot(alpha[stall_index] * 180/np.pi, cl[stall_index], 'ro', label="Stall angle of attack: {:.2f}".format(alpha[stall_index] * 180/np.pi))
        plt.plot(alpha0 * 180/np.pi, 0, 'go', label="Zero lift angle of attack: {:.2f}".format(alpha0 * 180/np.pi))
        plt.plot(alpha * 180/np.pi, m0 * alpha + cl0, 'k--', label="Linear regression")
        plt.xlabel("Angle of attack [deg]")
        plt.ylabel("Coefficient")
        plt.legend()
        plt.grid()
        plt.show()

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
    def __init__(self, radius:float=0.6, no_blades: int = 2, tip_speed: float = 100, cord: float = 0.1, twist:int=0, airfoil: str = "NACA23012")-> None:
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
        self.area = np.pi * (self.radius**2)
        self.no_blades = no_blades
        self.tip_speed = tip_speed
        self.omega = self.tip_speed / self.radius
        self.cord = cord
        self.twist = twist
        self.airfoil = airfoil
        self.rho = 1.225 # kg/m^3

        self.tip_loss = True

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

    def get_clalpha(self, plot=True):
        '''
        Returns the lift coefficient of the propeller
        '''
        v = 1.5e-5 # Dynamic viscosity of air [m^2/s]
        re = (self.tip_speed) * self.cord / v # Reynolds number

        # Run xfoil
        alpha, cl, cd, self.clalpha, alpha0, stall_index = xfoil(re=re, airfoil=self.airfoil, AOA_start=-7, AOA_end=20, AOA_step=0.1)
    
    def solidity(self):
        '''
        Returns the solidity of the propeller
        '''
        self.sigma = self.no_blades * self.cord / ( math.pi * self.radius)

    def BEMT(self,inflow:float=0, n_elements=10000, thrust:float=9*9.82, plot:bool=False, max_iter:int=100000):
        '''
        Calculates the thrust and power of the propeller using the Blade Element Momentum Theory
        twist: int, 0 = constant, 1 = linear, 2 = ideal
        n_elements: int, number of elements to divide the propeller into
        tip_loss: bool, if True, tip loss is accounted for
        '''
        self.solidity() # Get sigma
        self.get_clalpha() # Get clalpha

        # Calculate the CT
        Ct = thrust / (self.rho * self.tip_speed**2 * self.area)

        # Create arrays for the elements
        r = np.linspace(1/n_elements , 1-1/n_elements, n_elements)

        # Calculate the theta_tip
        theta_tip = 4 * Ct / (self.sigma * self.clalpha) + np.sqrt(Ct/2) # 3.84
        print(theta_tip)

        # Calculate the twist
        if self.twist == 0: # straight twist
            theta = np.ones(n_elements) * 6 * Ct / (self.sigma * self.clalpha) + 3/2 * np.sqrt(Ct/2) # 3.24
        elif self.twist == 1: # linear twist
            theta = r * (-theta_tip) + 2 * theta_tip # Figure 3.6
        elif self.twist == 2: # Ideal twist
            theta = theta_tip / r
        else:
            raise ValueError("Invalid twist value")
    
        # Calculate lambda
        if self.tip_loss:
            lam_init = np.ones(n_elements)
            for i in range(max_iter):
                lam = lam_init
                phi = lam / r

                f = self.no_blades / 2 * ((1 - r) / (r * phi))
                F = (2 / np.pi) * np.arccos(np.exp(-f))

                lam_init = ((self.sigma * self.clalpha) / (16 * F)) * (np.sqrt(1 + (32 * F) / (self.sigma * self.clalpha) * theta * r) - 1) # 3.126

                if np.abs(np.max(np.abs(lam_init) - np.abs(lam))) < 1e-8:
                    lam = lam_init
                    break
        else:
            lam = self.sigma * self.clalpha / (16) * (np.sqrt(1 + (32) / (self.sigma * self.clalpha) * theta * r) - 1) # 3.126

        self.lam = lam

        V_c = inflow / (self.omega * self.radius)

        self.CT = np.sum(4 * (V_c + self.lam)**2 * r * 1/n_elements) # 3.9
        self.CP = np.sum(4 * (V_c + self.lam)**3 * r * 1/n_elements) # 3.9

        self.thrust = self.CT * self.rho * self.tip_speed**2 * self.area
        self.power = self.CP * self.rho * self.tip_speed**3 * self.area


        if plot:
            plt.figure()
            plt.plot(r, F, label="Twist: {}".format(self.twist))
            plt.xlabel("Radius")
            plt.ylabel("F")
            plt.legend()
            plt.grid()

            plt.figure()
            plt.plot(r, self.lam, label="Twist: {}".format(self.twist))
            plt.xlabel("Radius")
            plt.ylabel("Lambda")
            plt.legend()
            plt.grid()

    def hover(self, force):
        '''
        Calculates the thrust and power of the propeller in hover with the given force
        force: float, force in N
        outputs:
        thrust: float, thrust of the propeller in N
        power: float, power of the propeller in W
        '''
        self.BEMT(thrust=force)

        return self.thrust, self.power
    
    def forward_flight(self, force, speed: float):
        '''
        Calculates the thrust and power of the propeller in forward flight
        '''
        self.BEMT(inflow=speed, thrust=force)

        return self.thrust, self.power


class Wing:
    def __init__(self, min_lift:float = 18, min_speed:float = 11, max_speed:float = 22, aspect_ratio:float = 10, airfoil:str="NACA4415",optimum_speed:float = 16):
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

    def reynolds(self, speed: float, c = 0.5): # 
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

    def stall_dim(self, plot:bool=False):
        '''
        Calculates the dimensions from the stall angle of attack and minimum speed
        '''
        self.xfoil_data()
        self.coefficient_lift()

        CL_max = self.CL[self.stall_index]
        S = 2 * self.min_lift / (self.min_speed**2 * self.rho * CL_max) # Wing area
        b = np.sqrt(self.AR * S) # Wing span
        c0 = 4 * S / (np.pi * b) # Root chord

        if plot: # plot wing shape
            angles = np.linspace(0, 2*np.pi, 100)
            x = b/2 * np.cos(angles)
            y = c0/2 * np.sin(angles)
            plt.figure()
            plt.plot(x, y)
            plt.title("Wing shape")
            plt.xlabel("Span [m]")
            plt.ylabel("Chord [m]")
            plt.axis('equal')
            plt.grid()

            # np.savetxt("Assignment3/stall_dim.txt", np.array([x, y]).T)

        return S, b, c0


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

            # np.savetxt(f"Assignment3/optimal_clcd{self.AR}.txt", np.array([self.alpha*180/np.pi, CLCD]).T)
            np.savetxt("Assignment3/optimal_clcdpoint.txt", np.array([self.alpha_optimal*180/np.pi, self.CL_optimal/self.CD_optimal]))

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
            x = self.b/2 * np.cos(angles)
            y = self.c0/2 * np.sin(angles)
            plt.figure()
            plt.plot(x, y)
            plt.title("Wing shape")
            plt.xlabel("Span [m]")
            plt.ylabel("Chord [m]")
            plt.axis('equal')
            plt.grid()

            np.savetxt("Assignment3/dimensions.txt", np.array([x, y]).T)

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

            np.savetxt("Assignment3/maxspeed_clcd.txt", np.array([self.alpha_max*180/np.pi, self.CL[self.index_max]/self.CD[self.index_max]]))

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

            np.savetxt("Assignment3/minspeed_clcd.txt", np.array([self.alpha_min*180/np.pi, self.CL[self.index_min]/self.CD[self.index_min]]))

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
    if False: # Test
        prop = Prop(twist=2)
        prop.BEMT(inflow=10, thrust = 1 ,plot=True)
        print(1)
        print("CT: {:.6f}".format(prop.CT))
        print("CP: {:.6f}".format(prop.CP))
        print("Thrust: {:.6f}".format(prop.thrust))
        print("Power: {:.6f}".format(prop.power))
        plt.show()

    if False: # Test
        plane = Plane()
        plane.fixed_aoa(plot=True)

        plt.show()

    if False: # Question 1
        kappa = [0.6, 0.8, 1, 1.2, 1.4]
        R = np.linspace(0.2, 3, 100)
        T = 177/2 # N
        rho = 1.225 # kg/m^3
        N_blades = 2
        V_tip = 100 # m/s
        c = 0.05 # m
        Cd = 0.05

        plt.figure()
        for k in kappa:
            P = k * T**(3/2) / (np.sqrt(2 * np.pi * rho) * R) + (rho * N_blades * c * Cd * V_tip**3) / 8
            plt.plot(R, P, label=f"kappa = {k}")        
        plt.plot(R, P)
        plt.xlabel("Radius [m]")
        plt.ylabel("Power [W]")
        plt.grid()
        plt.legend()
        plt.show()



    if False: # Question 3
        props = [Prop(twist=0), Prop(twist=1), Prop(twist=2)]

        Thrust = np.linspace(9*9.8, 12*9.8, 100)
        thrust = np.zeros(len(Thrust))
        power = np.zeros(len(Thrust))
        plt.figure()
        for i, prop in enumerate(props):
            for j, force in enumerate(Thrust):
                thrust[j], power[j] = prop.hover(force)
            plt.plot(Thrust, power, label=f"Twist: {i}")
        plt.xlabel("Thrust [N]")
        plt.ylabel("Power [W]")
        plt.legend()
        plt.grid()
        plt.show()

    if True: # Question 4
        if False: # Hurtige måde
            wing = Wing()
            print(wing.stall_dim(plot=True))

        else:
            wing = Wing(aspect_ratio=10)
            wing.xfoil_data(plot=False)
            wing.optimal_ClCd(plot=True)
            wing.maxspeed_ClCd(plot=True)
            wing.minspeed_ClCd(plot=True)
            print("optimal angle of attack: {:.6f}".format(wing.alpha_optimal * 180/np.pi))
            print("optimal Cl/Cd: {:.6f}".format(wing.CL_optimal / wing.CD_optimal))
            print("optimal Cl: {:.6f}".format(wing.CL_optimal))
            print("optimal Cd: {:.6f}".format(wing.CD_optimal))
            print("Max speed angle of attack: {:.2f}".format(wing.alpha_max * 180/np.pi))
            print("Max speed Cl: {:.6f}".format(wing.CL_max))
            print("Max speed Cd: {:.6f}".format(wing.CD_max))
            wing.dimensions(plot=True)
            print(wing.S)
            print(wing.b)
            print(wing.c0)

        plt.show()

        print("stop")

    if False: # Question 5
        wing = Wing(min_lift=18, min_speed=11, max_speed=20, aspect_ratio=10)
        wing.CLCD_speed(plot=True)

        plt.show()