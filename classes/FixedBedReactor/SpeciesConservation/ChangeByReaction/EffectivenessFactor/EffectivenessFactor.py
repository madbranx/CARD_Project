# Class to calculate the effectiveness factor based on the thiele modulus for spherical cat particles
# Till Kasselmann 16.01.2025

import numpy as np
import casadi as casADi

class EffectivenessFactor:
    def __init__(self, log, RSQ, GCF):
        self.log = log
        self.RSQ = RSQ
        self.GCF = GCF
        self.R = RSQ.getParameterValue("R")
        self.pi = np.pi

    def calc(self, T, w_i):

        thiele = self.__calc_thiele(T, w_i)
        effectiveness_factor = 3/thiele * ( 1/casADi.tanh(thiele) - 1/thiele )
        return effectiveness_factor

    def __calc_thiele(self, T, w_i):

        # TODO calc from w_i array
        concentration_CO2 = 1
        stochiometry_CO2 = self.RSQ.getStoichCoeffs()[1]
        diameter_particle = self.RSQ.getParameterValue("diameter_particle")

        # TODO calc from w_i array
        [p_CH4, p_H2O, p_CO2, p_H2 ] = [1, 2, 3, 4]

        eff_diff_coff = self.__calc_eff_diff_coff(T)

        #TODO
        # reaction rate wie impelemtieren??
        # fix call
        reaction_rate = self.RSQ.getReactionRate()
        reaction_rate = reaction_rate(T, p_CH4, p_H2O, p_CO2, p_H2)

        thiele = diameter_particle / 2 * ((stochiometry_CO2 * reaction_rate) / (eff_diff_coff * concentration_CO2)) ** (0.5)
        return thiele

    def __calc_eff_diff_coff(self, T):
        tortuosity_particle = self.RSQ.getParameterValue("tortuosity_particle")
        porosity_particle = self.RSQ.getParameterValue("porosity_particle")


        # Checked
        knudsen_diff_coff = self.__calc_knudsen_diff_coff(T)

        #TODO TEMP
        # molar diff coff, wie implementieren??
        molar_diff_coff = 7.1497e-5


        eff_diff_coff = ( tortuosity_particle**2 / porosity_particle * (1/molar_diff_coff + 1/knudsen_diff_coff) )**(-1)
        return eff_diff_coff

    def __calc_knudsen_diff_coff(self, T):
        # Checked
        diameter_pore = self.RSQ.getParameterValue("diameter_pore")
        CO2 = self.RSQ.getComponent("CO2")
        molar_mass_i = CO2.get_molar_mass()

        knudsen_diff_coff = diameter_pore / 3 * ( (8 * self.R * T) / (self.pi * molar_mass_i) )**(0.5)
        return knudsen_diff_coff