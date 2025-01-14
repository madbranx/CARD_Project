# Class to calculate the kinetics
# Till Kasselmann 13.01.2025

import casadi as casADi

class ReactionRate:
    def __init__(self):
        self.T_ref = 550        # K
        self.k0_ref = 3.46e-4   # mol/(bar s g_cat)
        self.E_A = 77.5e3       # J/mol
        self.R = 8.314          # J/(mol*K)
        self.AdsorptionConstants = {
            'OH': {'K_ref': 0.5, 'dH': 22.4e3},     # bar^(-0.5), J/mol
            'H2': {'K_ref': 0.44, 'dH': -6.2e3},    # bar^(-0.5), J/mol
            'mix': {'K_ref': 0.88, 'dH': -10e3},    # bar^(-0.5), J/mol
        }

    # methode to calculate rate coefficient k
    def __calc_k (self, T):
        # T in K
        k = self.k0_ref * casADi.exp(self.E_A / self.R *( 1/self.T_ref - 1/T))
        return k

    # methode to calculate adsorption constant K_x for OH, H2 and mix
    def __calc_K_x(self, species, T):
        # T in K
        # Species string input for AdsorptionConstants dictionary
        if species not in self.AdsorptionConstants:
            raise ValueError(f"Unknown species in adsorption constant K_x: {species}")

        K_x_ref = self.AdsorptionConstants[species]['K_ref']
        dH_x = self.AdsorptionConstants[species]['dH']

        K_x = K_x_ref * casADi.exp(dH_x/self.R *( 1/self.T_ref - 1/T))
        return K_x

    # methode to calculate equilibrium constant K_eq
    def __K_eq(self, T):
        # T in K
        K_eq = 137 * T**(-3.998) * casADi.exp(158.7e3 / (self.R * T))
        return K_eq

    def rate_equation(self, T, p_CH4, p_H2O, p_CO2, p_H2, rho_cat):
        # T in K
        # p_i in bar
        # rho_cat in g_cat / m^3_cat
        k = self.__calc_k(T)
        K_eq = self.__K_eq(T)

        r = (k * (p_CO2**0.5) * (p_H2**0.5) * (1 - (p_CH4 * (p_H2O**2)) / (K_eq * p_CO2 * (p_H2**4))) /
             (1 + self.__calc_K_x('OH', T) * p_H2O / (p_H2 ** 0.5) + self.__calc_K_x('H2', T) * (p_H2 ** 0.5) + self.__calc_K_x('mix', T) * (p_CO2 ** 0.5)))

        # Conversion from mol/(g_cat*s) into mol/(m^3_cat*s)
        r = r * rho_cat * 1000
        return r

