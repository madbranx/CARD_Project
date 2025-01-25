from classes.Parameters.Component import Component, MaterialProperty

import numpy as np

class Parameters:
    def __init__(self):
        # add Constants to RSQ
        self.R =  8.314
        self.pi = np.pi

        # add Parameters to RSQ
        self.reactorLength=2.5
        self.reactorDiameter=0.02

        self.cat_diameter=0.002
        self.cat_tortuosity=2
        self.cat_porosity=0.6

        self.diameter_pore=10e-9

        #self.eps = self.__calculate_void_fraction()
        self.eps = 0.4 # given by task

        self.T_in=300
        self.u_in=1
        self.w_i_in=[0, 0, 0.2, 0.8] # CH4, H20, CO2, H2
        self.p_in=5e5

        # Parameters given by task for 1D calculation
        self.T_wall=550
        self.lambda_radial=200

        # Components (incl. cat)
        self.components = None
        self.Cat = None
        self.__initializeComponents()

    def __initializeComponents(self):
        # add Components and their Properties to RSQ
        CH4 = Component("CH4") # SOURCE: NIST-Database with polynomial fit at 5 bar
        CH4.add_property(Component.DENSITY, 0.348)
        CH4.add_property(Component.HEAT_CAPACITY, [25.455571, 2.549214e-02, 3.318571e-05], MaterialProperty.POLYNOMIAL)
        CH4.add_property(Component.THERMAL_CONDUCTIVITY, [2.054738e-03, 7.050143e-05, 1.239619e-07], MaterialProperty.POLYNOMIAL)
        CH4.add_property(Component.COLLISION_AREA, 0.46e-18)
        CH4.add_property(Component.DIFFUSION_VOLUME, 25.14)
        CH4.add_property(Component.DYNAMIC_VISCOSITY, [5.223333e-07, 3.968143e-08, -1.364762e-11], MaterialProperty.POLYNOMIAL)
        CH4.add_property(Component.MOLECULAR_WEIGHT, 0.01604)

        H20 = Component("H2O") # SOURCE: NIST-Database with polynomial fit at 5 bar
        H20.add_property(Component.DENSITY, 0.322)
        H20.add_property(Component.HEAT_CAPACITY, [43.064143, -2.033657e-02, 1.881429e-05], MaterialProperty.POLYNOMIAL)
        H20.add_property(Component.THERMAL_CONDUCTIVITY, [-2.047571e-03, 5.770379e-05, 4.061786e-08], MaterialProperty.POLYNOMIAL)
        H20.add_property(Component.COLLISION_AREA, 0.46e-18)
        H20.add_property(Component.DIFFUSION_VOLUME, 13.1)
        H20.add_property(Component.DYNAMIC_VISCOSITY, [-4.954286e-06, 4.593789e-08, -3.341071e-12], MaterialProperty.POLYNOMIAL)
        H20.add_property(Component.MOLECULAR_WEIGHT, 0.01802)

        CO2 = Component("CO2") # SOURCE: NIST-Database with polynomial fit at 5 bar
        CO2.add_property(Component.DENSITY, 0.725)
        CO2.add_property(Component.HEAT_CAPACITY, [27.000304, 4.421494e-02, -1.689345e-05], MaterialProperty.POLYNOMIAL)
        CO2.add_property(Component.THERMAL_CONDUCTIVITY, [-9.104310e-03, 8.906750e-05, -8.951190e-09], MaterialProperty.POLYNOMIAL)
        CO2.add_property(Component.COLLISION_AREA, 0.52e-18)
        CO2.add_property(Component.DIFFUSION_VOLUME, 26.9)
        CO2.add_property(Component.DYNAMIC_VISCOSITY, [-6.721429e-08, 5.467619e-08, -1.347381e-11], MaterialProperty.POLYNOMIAL)
        CO2.add_property(Component.MOLECULAR_WEIGHT, 0.04401)

        H2 = Component("H2") # SOURCE: NIST-Database with polynomial fit at 5 bar
        H2.add_property(Component.DENSITY, 0.041)
        H2.add_property(Component.HEAT_CAPACITY, [28.838423, 1.066071e-04, 1.225595e-06], MaterialProperty.POLYNOMIAL)
        H2.add_property(Component.THERMAL_CONDUCTIVITY, [6.286452e-02, 4.312905e-04, -3.483333e-08], MaterialProperty.POLYNOMIAL)
        H2.add_property(Component.COLLISION_AREA, 0.27e-18)
        H2.add_property(Component.DIFFUSION_VOLUME, 6.12)
        H2.add_property(Component.DYNAMIC_VISCOSITY, [2.801093e-06, 2.172896e-08, -3.832798e-12], MaterialProperty.POLYNOMIAL)
        H2.add_property(Component.MOLECULAR_WEIGHT, 0.002016)

        self.components = [CH4, H20, CO2, H2]

        # Catalyst Properties
        cat = Component("cat")
        cat.add_property(Component.DENSITY, 2355.2)
        cat.add_property(Component.HEAT_CAPACITY, 1107)
        cat.add_property(Component.THERMAL_CONDUCTIVITY, 3.6)

        self.cat = cat

    def __calculate_void_fraction(self):
        ratio = self.cat_diameter / self.reactorDiameter

        if ratio <= 0.5:
            epsilon = 0.4 + 0.05 * ratio + 0.412 * (ratio ** 2)
        elif 0.5 < ratio <= 0.536:
            epsilon = 0.528 + 2.464 * (ratio - 0.5)
        else:  # ratio >= 0.536
            epsilon = 1 - 0.667 * (ratio ** 3) * (2 * ratio - 1) ** -0.5

        return epsilon
