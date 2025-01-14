# Class to calculate axial mass flow
# Till Kasselmann 14.01.2025

class AxialMassFlow:

    def __init__(self):
        pass

    def calc(self,axial_flowfield, mass_fraction, fluid_density):
        # axial_flowfield in m/s
        # fluid_density in kg/m^3
        axial_mass_flux = - axial_flowfield * mass_fraction * fluid_density
        return axial_mass_flux
