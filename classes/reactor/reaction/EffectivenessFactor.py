# Class to calculate the effectiveness factor based on the thiele modulus for spherical cat particles
# Till Kasselmann 16.01.2025
from bdb import effective

import numpy as np
import casadi as casADi

class EffectivenessFactor:
    def __init__(self):
        self.pi = np.pi
        self.R = 8.314

    def calc_effectiveness_factor(self, T, molar_mass_i, stochiometry_i, concentration_i, diameter_pore, tortuosity_particle, porosity_particle, diameter_particle):
        thiele = self.__calc_thiele(T, molar_mass_i, stochiometry_i, concentration_i, diameter_pore, tortuosity_particle, porosity_particle, diameter_particle)
        effectiveness_factor = 3/thiele * ( 1/casADi.tanh(thiele) - 1/thiele )
        return effectiveness_factor

    def __calc_thiele(self, T, molar_mass_i, stochiometry_i, concentration_i, diameter_pore, tortuosity_particle, porosity_particle, diameter_particle):
        eff_diff_coff = self.__calc_eff_diff_coff(T, molar_mass_i, diameter_pore, tortuosity_particle, porosity_particle)
        # reaction rate wie impelemtieren??
        reaction_rate = 1   #TEMP
        thiele = diameter_particle / 2 * ( (stochiometry_i * reaction_rate) / (eff_diff_coff * concentration_i) )**(0.5)
        return thiele

    def __calc_eff_diff_coff(self, T, molar_mass_i, diameter_pore, tortuosity_particle, porosity_particle):
        knudsen_diff_coff = self.__calc_knudsen_diff_coff(T, molar_mass_i, diameter_pore)
        # molar diff coff, wie implementieren??
        molar_diff_coff = 1 # TEMP
        eff_diff_coff = ( tortuosity_particle**2 / porosity_particle * (1/molar_diff_coff + 1/knudsen_diff_coff) )**(-1)
        return eff_diff_coff

    def __calc_knudsen_diff_coff(self, T, molar_mass_i, diameter_pore):
        knudsen_diff_coff = diameter_pore / 3 * ( (8 * self.R * T) / (self.pi * molar_mass_i) )**(0.5)
        return knudsen_diff_coff