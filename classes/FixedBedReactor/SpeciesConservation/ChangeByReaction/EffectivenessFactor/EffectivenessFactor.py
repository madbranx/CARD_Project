# Class to calculate the effectiveness factor based on the thiele modulus for spherical cat particles
# Till Kasselmann 16.01.2025

import numpy as np
import casadi as casADi

class EffectivenessFactor:
    def __init__(self, log, RSQ, GCF):
        self.log = log
        self.RSQ = RSQ
        self.GCF = GCF

    def calc(self, w_i, T, p):

        thiele = self.__calc_thiele(T, w_i, p)
        effectiveness_factor = 3/thiele * ( 1/casADi.tanh(thiele) - 1/thiele)
        return effectiveness_factor

    def __calc_thiele(self, T, w_i, p):
        c_i = self.GCF.concentrations(w_i, T, p)
        concentration_CO2 = c_i[2]

        diameter_particle = self.RSQ.getParameterValue("cat_diameter")

        # get the relevant reaction for effectiveness factor
        reaction_1 = self.RSQ.getReactions()[0]
        stoichiometry_coefficients = reaction_1.getStoichiometryCoefficients()
        stoichiometry_CO2 = stoichiometry_coefficients[2]
        reaction_rate = reaction_1.getReactionRate(w_i, T, p)

        # calculate the effective diffusion coefficient

        eff_diff_coff = self.__calc_eff_diff_coff(T)

        thiele = diameter_particle / 2 * ((stoichiometry_CO2 * reaction_rate) / (eff_diff_coff * concentration_CO2)) ** (0.5)
        return thiele

    def __calc_eff_diff_coff(self, T):
        tortuosity_particle = self.RSQ.getParameterValue("cat_tortuosity")
        porosity_particle = self.RSQ.getParameterValue("cat_porosity")

        knudsen_diff_coff = self.__calc_knudsen_diff_coff(T)

        #TODO
        # molar diff coff, wie implementieren??
        molar_diff_coff = 7.1497e-5

        eff_diff_coff = ( tortuosity_particle**2 / porosity_particle * (1/molar_diff_coff + 1/knudsen_diff_coff) )**(-1)
        return eff_diff_coff

    def __calc_knudsen_diff_coff(self, T):
        diameter_pore = self.RSQ.getParameterValue("diameter_pore")
        pi = np.pi
        R = self.RSQ.getParameterValue("R")

        CO2 = self.RSQ.getComponent("CO2")
        molar_mass_i = CO2.get_molecular_weight()

        knudsen_diff_coff = diameter_pore / 3 * ( (8 * R * T) / (pi * molar_mass_i) ) ** 0.5
        return knudsen_diff_coff