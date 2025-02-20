from classes.Parameters.Component import Component, MaterialProperty

import numpy as np

class Parameters:
    def __init__(self):

        # Use of constant values @500 K
        self.temperature_dependent_matProps = True

        self.R =  8.314  # J/(kgK)
        self.pi = np.pi
        self.boltzmann = 1.380649e-23    # J/K
        self.radiation_blackBody = 5.67e-8 # W/(m^2 K^4)

        # add Parameters
        self.reactorLength = 2.5  # m
        self.reactorDiameter = 0.02  # m

        self.cat_diameter = 0.002  # m
        self.cat_tortuosity = 2  # -
        self.cat_porosity = 0.6  # -
        self.cat_emissionCoefficient = 0.9  # -

        self.reactor_thermalConductivity = 20   # W/(m K)    # Bremer
        self.reactor_wallThickness = 0.002       # m   # Bremer


        self.diameter_pore=10e-9

        self.eps = self.calculate_void_fraction() # = 0.409 for 0.002 cat diameter
        #self.eps = 0.4 # given by task

        self.T_in = 300  # K
        self.u_in= 1  # m/s
        # molar inlet ratio H2:CO2 = 4:1
        self.w_i_in = [0, 0, 0.8452, 0.1548] # kg/kg CH4, H20, CO2, H2
        self.p_in = 5e5 # Pa

        # Parameters given by task for 1D calculation
        self.T_wall = 550  # K
        self.lambda_radial = 200  # W/(m^2K)

        # Components (incl. cat)
        self.components = None
        self.Cat = None

        self.__initializeComponents()

    def __initializeComponents(self):
        ## ADD COMPONENTS AND THEIR PROPERTIES
        # Density                       kg/m^3
        # Heat Capacity                 J/(kg K)
        # Thermal Conductivity          W/(m K)
        # Dynamic Viscosity             Pa s
        # Molecular Weight              kg/mol
        # Diffusion Volume              -

    # CH4
        CH4 = Component("CH4") # SOURCE: NIST-Database with polynomial fit at 5 bar
        CH4.add_property(Component.DENSITY, 0.348)

        if self.temperature_dependent_matProps is True:
            CH4.add_property(Component.HEAT_CAPACITY, [1587.0057,1.5893,0.00207], MaterialProperty.POLYNOMIAL)
            CH4.add_property(Component.THERMAL_CONDUCTIVITY, [2.054738e-03, 7.050143e-05, 1.239619e-07], MaterialProperty.POLYNOMIAL)
            CH4.add_property(Component.DYNAMIC_VISCOSITY, [5.223333e-07, 3.968143e-08, -1.364762e-11], MaterialProperty.POLYNOMIAL)
        else:
            CH4.add_property(Component.HEAT_CAPACITY, 2899.1557)
            CH4.add_property(Component.THERMAL_CONDUCTIVITY, 0.068295928)
            CH4.add_property(Component.DYNAMIC_VISCOSITY, 1.69511433e-05)

        CH4.add_property(Component.MOLECULAR_WEIGHT, 0.01604)
        CH4.add_property(Component.DIFFUSION_VOLUME, 25.14)

    # H2O
        H2O = Component("H2O") # SOURCE: NIST-Database with polynomial fit at 5 bar
        H2O.add_property(Component.DENSITY, 0.322)

        if self.temperature_dependent_matProps is True:
            H2O.add_property(Component.HEAT_CAPACITY, [2389.797, -1.129, 0.001044], MaterialProperty.POLYNOMIAL)
            H2O.add_property(Component.THERMAL_CONDUCTIVITY, [-2.047571e-03, 5.770379e-05, 4.061786e-08], MaterialProperty.POLYNOMIAL)
            H2O.add_property(Component.DYNAMIC_VISCOSITY, [-4.954286e-06, 4.593789e-08, -3.341071e-12], MaterialProperty.POLYNOMIAL)
        else:
            H2O.add_property(Component.HEAT_CAPACITY, 2086.297)
            H2O.add_property(Component.THERMAL_CONDUCTIVITY, 0.036958789)
            H2O.add_property(Component.DYNAMIC_VISCOSITY, 1.717939125e-05)

        H2O.add_property(Component.MOLECULAR_WEIGHT, 0.01802)
        H2O.add_property(Component.DIFFUSION_VOLUME, 13.1)

    # CO2
        CO2 = Component("CO2") # SOURCE: NIST-Database with polynomial fit at 5 bar
        CO2.add_property(Component.DENSITY, 0.725)

        if self.temperature_dependent_matProps is True:
            CO2.add_property(Component.HEAT_CAPACITY, [613.504,1.005,-0.000384], MaterialProperty.POLYNOMIAL)
            CO2.add_property(Component.THERMAL_CONDUCTIVITY, [-9.104310e-03, 8.906750e-05, -8.951190e-09], MaterialProperty.POLYNOMIAL)
            CO2.add_property(Component.DYNAMIC_VISCOSITY, [-6.721429e-08, 5.467619e-08, -1.347381e-11], MaterialProperty.POLYNOMIAL)
        else:
            CO2.add_property(Component.HEAT_CAPACITY, 1020.004)
            CO2.add_property(Component.THERMAL_CONDUCTIVITY, 0.0331916425)
            CO2.add_property(Component.DYNAMIC_VISCOSITY, 2.390242821e-05)

        CO2.add_property(Component.MOLECULAR_WEIGHT, 0.04401)
        CO2.add_property(Component.DIFFUSION_VOLUME, 26.9)

    # H2
        H2 = Component("H2") # SOURCE: NIST-Database with polynomial fit at 5 bar
        H2.add_property(Component.DENSITY, 0.041)

        if self.temperature_dependent_matProps is True:
            H2.add_property(Component.HEAT_CAPACITY, [14304.773, 0.05288, 0.000608], MaterialProperty.POLYNOMIAL)
            H2.add_property(Component.THERMAL_CONDUCTIVITY, [6.286452e-02, 4.312905e-04, -3.483333e-08], MaterialProperty.POLYNOMIAL)
            H2.add_property(Component.DYNAMIC_VISCOSITY, [2.801093e-06, 2.172896e-08, -3.832798e-12], MaterialProperty.POLYNOMIAL)
        else:
            H2.add_property(Component.HEAT_CAPACITY, 14483.213)
            H2.add_property(Component.THERMAL_CONDUCTIVITY, 0.2698014375)
            H2.add_property(Component.DYNAMIC_VISCOSITY, 1.27073735)

        H2.add_property(Component.MOLECULAR_WEIGHT, 0.002016)
        H2.add_property(Component.DIFFUSION_VOLUME, 6.12)

        # check_temperature = 500
        # self.checkPropertyValues(CH4, check_temperature)
        # self.checkPropertyValues(H2O, check_temperature)
        # self.checkPropertyValues(CO2, check_temperature)
        # self.checkPropertyValues(H2, check_temperature)

        self.components = [CH4, H2O, CO2, H2]


    # CAT
        # Catalyst Properties
        cat = Component("cat")
        cat.add_property(Component.DENSITY, 2355.2)
        cat.add_property(Component.HEAT_CAPACITY, 1107)
        cat.add_property(Component.THERMAL_CONDUCTIVITY, 3.6)

        self.cat = cat

    def calculate_void_fraction(self):
        ratio = self.cat_diameter / self.reactorDiameter

        if ratio <= 0.5:
            epsilon = 0.4 + 0.05 * ratio + 0.412 * (ratio ** 2)
        elif 0.5 < ratio <= 0.536:
            epsilon = 0.528 + 2.464 * (ratio - 0.5)
        else:  # ratio >= 0.536
            epsilon = 1 - 0.667 * (ratio ** 3) * (2 * ratio - 1) ** -0.5

        return epsilon


    def checkPropertyValues(self, comp, temperature):
        print("Material Properties of " + comp.get_Name() + " @T=" + str(temperature) + " K")

        print("Density", comp.get_property(Component.DENSITY, temperature), " kg/m^3")
        print("Heat Capacity", comp.get_property(Component.HEAT_CAPACITY, temperature), " J/(kg K)")
        print("Thermal Conductivity", comp.get_property(Component.THERMAL_CONDUCTIVITY, temperature), " W/(m K)")
        print("Dynamic Viscosity", comp.get_property(Component.DYNAMIC_VISCOSITY, temperature), " Pa s")

        print("Molecular Weight", comp.get_property(Component.MOLECULAR_WEIGHT, temperature), " kg/mol")
        print("Diffusion Volume", comp.get_property(Component.DIFFUSION_VOLUME, temperature), "\n")
