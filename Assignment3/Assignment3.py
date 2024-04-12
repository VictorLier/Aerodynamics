import numpy as np
import math

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
    def __init__(self, diameter: float, no_blades: int, tip_speed: float):
        '''
        self.area = area
        self.no_blades = no_blades
        self.tip_speed = tip_speed
        '''
        self.diameter = diameter
        self.area = np.pi * (self.diameter / 2)**2
        self.no_blades = no_blades
        self.tip_speed = tip_speed
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
        battery = Battery(mass = 4)
        print("Battery mass =",battery.mass)
        print("Battery energy =",battery.energy)
        print("Battery power =",battery.power)

        prop = Prop(diameter = 0.3048, no_blades = 2, tip_speed = 30)
        print("Prop ideal_thrust =",prop.ideal_thrust())
        print("Prop ideal_power =",prop.ideal_power())

        plane = Plane(mass = 18, battery = battery, hover_prop = prop, tilt_prop = prop)
        print("Plane payload =",plane.payload())
        print("Plane ideal_hover_time =",plane.ideal_hover_time())
        print("Plane required_hover_thrust =",plane.required_hover_thrust())
        print("Prop ideal_thrust =",prop.ideal_thrust())