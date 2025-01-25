from V2.classes.Parameters.Component import Component
from V2.classes.Properties.Kinetics import Kinetics

import casadi as CasADi

class SpeciesConservation(Kinetics):
    def __init__(self):
        super().__init__()

    ## AXIAL MASS FLOW
    def axialMassFlow(self, T, w_i, w_i_in, u, p, comp):
        # Only convection
        rho_fl = self.rho_fl(w_i, T, p)

        j_i_ax = -u*rho_fl*(w_i[comp]-w_i_in[comp])
        return j_i_ax

    ## RADIAL MASS FLOW #TODO
    def radialMassFlow(self, T, w_i, w_i_in, u, p, comp):
        pass

    ## CHANGE BY REACTION
    def changeByReaction(self, T, w_i, p, comp):
        Mw_i = self.getMolarWeights()
        eff_factor = self.effFactor(w_i, T, p)
        reaction_rate = self.rate_equation(w_i, T, p)

        return (1 - self.eps) * Mw_i[comp] * self.nu[comp] * eff_factor * reaction_rate


    def effFactor(self, w_i, T, p):
        # effectiveness factor for CO2 used
        # index of CO2 = 2
        comp = 2
        thiele = self.__calc_thiele(T, w_i, p, comp)
        effectiveness_factor = 3 / thiele * (1 / CasADi.tanh(thiele) - 1 / thiele)
        return effectiveness_factor

    def __calc_thiele(self, T, w_i, p, comp):
        c_i = self.concentrations(w_i, T, p)
        concentration_CO2 = c_i[2]
        stoichiometry_CO2 = self.nu[2]
        diameter_particle = self.cat_diameter
        reaction_rate = self.rate_equation(w_i, T, p)

        # calculate the effective diffusion coefficient
        eff_diff_coff = self.calc_eff_diff_coff(T, w_i, p, comp)

        thiele = diameter_particle / 2 * ((stoichiometry_CO2 * reaction_rate) / (eff_diff_coff * concentration_CO2)) ** 0.5
        return thiele

    def calc_eff_diff_coff(self, T, w_i, p, comp):
        tortuosity_particle = self.cat_tortuosity
        porosity_particle = self.cat_porosity
        molar_mass = self.getMolarWeights()[comp]

        knudsen_diff_coff = self.Knudsen_diff_coff(T, molar_mass)

        molar_diff_coff = self.MixtureAveragedDiffusionCoefficient(w_i, T, p, comp)

        eff_diff_coff = 1 / ((tortuosity_particle ** 2 / porosity_particle) * (
                    1 / molar_diff_coff + 1 / knudsen_diff_coff))
        return eff_diff_coff

    def Knudsen_diff_coff(self, T, molarMass):

        knudsen_diff_coff = self.diameter_pore / 3 * ((8 * self.R * T) / (self.pi * molarMass)) ** 0.5
        return knudsen_diff_coff

    def MixtureAveragedDiffusionCoefficient(self, w_i, T, p, comp): # i = 0,1,.... depending on the component index
        # get diffusion volumes, molar masses and densities for all components
        diff_volumes = []
        molar_weights = []
        densities = []

        for component in self.components:
            diff_volumes.append(component.get_property(Component.DIFFUSION_VOLUME))
            molar_weights.append(component.get_property(Component.MOLECULAR_WEIGHT))
            #densities.append(component.get_property(Component.DENSITY, T)) # Using material properties

        densities = self.rho_comp(w_i, T, p) # Using ideal Gas Law

        # get fluid mixture values for the density and molar mass
        rho_fl = self.rho_fl(w_i, T, p)
        Mw_fl = self.massFraction_weighted_average(w_i, Component.MOLECULAR_WEIGHT)

        # calculating sum of rho_j/ Mw_j * binary diffusion coefficients_i,j
        summ = 0
        for j in range(len(self.components)):
            if j != comp:
                fuller_ij = self.Fuller(T, p, molar_weights[comp], molar_weights[j], diff_volumes[comp], diff_volumes[j])
                summ += densities[j] / (fuller_ij * molar_weights[j])

        avgDiffusionCoefficient = ((rho_fl/Mw_fl) - (densities[comp] / molar_weights[comp])) / summ
        return avgDiffusionCoefficient

    def Fuller(self, T, p, molar_mass_i, molar_mass_j, diff_volume_i, diff_volume_j):
        # Conversion from SI units to equation units
        p = p * 1e-5                           # Pa -> bar
        molar_mass_i = molar_mass_i * 1e3      # kg/mol -> g/mol
        molar_mass_j = molar_mass_j * 1e3      # kg/mol -> g/mol

        # Diffusion coefficient of i in j after Fuller in m^2/s
        DiffCoff = (1e-4 * (0.00143 * T ** 1.75 * ((1/molar_mass_i) + (1/molar_mass_j)) ** 0.5) /
                    (p * 2 ** 0.5 * (diff_volume_i ** (1 / 3) + diff_volume_j ** (1 / 3)) ** 2))

        return DiffCoff
