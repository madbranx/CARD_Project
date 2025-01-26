import casadi as CasADi

from classes.Parameters.Component import Component
from classes.Properties.FluidProperties import FluidProperties


class Kinetics(FluidProperties):
    def __init__(self):
        super().__init__()

        self.nu = [-1, -2, 1, 4]

        # Kinetic specific Variables
        self.T_ref = 550        # K
        self.k0_ref = 3.46e-4   # mol/(bar s g_cat)
        self.E_A = 77.5e3       # J/mol
        self.AdsorptionConstants = {
            'OH': {'K_ref': 0.5, 'dH': 22.4e3},     # bar^(-0.5), J/mol
            'H2': {'K_ref': 0.44, 'dH': -6.2e3},    # bar^(-0.5), J/mol
            'mix': {'K_ref': 0.88, 'dH': -10e3},    # bar^(-0.5), J/mol
        }
        self.reactionEnthalpy = -164900

    def rate_equation(self, w_i, T, p):
        # p_i in bar
        # [p_CH4, p_H2O, p_CO2, p_H2] = self.partial_pressures(w_i, p) * 1e-5
        p_i = self.partial_pressures(w_i, p) * 1e-5
        p_CH4 = p_i[0]
        p_H2O = p_i[1]
        p_CO2 = p_i[2]
        p_H2 = p_i[3]

        k = self.__calc_k(T)
        K_eq = self.__K_eq(T)

        r = (k * (p_CO2 ** 0.5) * (p_H2 ** 0.5) * (1 - (p_CH4 * (p_H2O ** 2)) / (K_eq * p_CO2 * (p_H2 ** 4))) /
             (1 + self.__calc_K_x('OH', T) * p_H2O / (p_H2 ** 0.5) + self.__calc_K_x('H2', T) * (
                         p_H2 ** 0.5) + self.__calc_K_x('mix', T) * (p_CO2 ** 0.5))**2)

        # rho_cat from [kg_cat / m^3_cat] into [g_cat / m^3_cat]
        rho_cat = self.cat.get_property(Component.DENSITY, T) * 1000

        # Conversion of r from [mol/(g_cat*s)] into [mol/(m^3_cat*s)]
        r *= rho_cat

        return r

    # methode to calculate rate coefficient k
    def __calc_k(self, T):
        # T in K
        k = self.k0_ref * CasADi.exp(self.E_A / self.R * (1 / self.T_ref - 1 / T))
        return k

    # methode to calculate adsorption constant K_x for OH, H2 and mix
    def __calc_K_x(self, species, T):
        # T in K
        K_x_ref = self.AdsorptionConstants[species]['K_ref']
        dH_x = self.AdsorptionConstants[species]['dH']

        K_x = K_x_ref * CasADi.exp(dH_x / self.R * (1 / self.T_ref - 1 / T))
        return K_x

    # methode to calculate equilibrium constant K_eq
    def __K_eq(self, T):
        # T in K
        K_eq = 137 * T ** (-3.998) * CasADi.exp(158.7e3 / (self.R * T))
        return K_eq
