from classes.Parameters.Parameters import Parameters
from classes.Parameters.Component import Component

import casadi as CasADi

"""
The FluidProperties class is used to define methods for calculation of parameters of the fluid phase.
This includes dimensionless quantities (Re, Pr, Pe), the thermal conductivity with viscosity mixing rule and general
properties as density and mass fraction weighted averages.
"""


class FluidProperties(Parameters):

    def __init__(self):
        super().__init__()

    def rho_comp(self, w_i, T, p):
        # Calculate density of components
        return self.concentrations(w_i, T, p)*self.getMolarWeights()

    def rho_fl(self, w_i, T, p):
        # Calculate density of fluid phase
        Mw_fl = self.massFraction_weighted_average(w_i, Component.MOLECULAR_WEIGHT)
        rho = self.__c_fl(T, p) * Mw_fl
        return rho

    def __c_fl(self, T, p):
        # Calculate mean fluid concentration
        cf = p / (self.R * T)
        return cf

    def getMolarWeights(self):
        Mw_i = []
        for component in self.components:
            Mw_i.append(component.get_property(Component.MOLECULAR_WEIGHT))
        return Mw_i

    def massFraction_weighted_average(self, w_i, material_property, T=None):
        # Methode to calculate mass fraction weighted averages of fluid properties
        # property = Component.DENSITY, ...
        property_values = []
        for component in self.components:
            property_values.append(component.get_property(material_property, T))
        summ = 0
        for i, property_value in enumerate(property_values):
            if material_property == Component.MOLECULAR_WEIGHT:
                summ += w_i[i] / property_value
            else:
                summ += w_i[i] * property_value
        if material_property == Component.MOLECULAR_WEIGHT:
            return 1 / summ
        else:
            return summ


    def moleFractions(self, w_i):
        # Methode to convert mass fractions into mole fractions
        Mw_i = self.getMolarWeights()
        moles = w_i / Mw_i
        tot_moles = CasADi.sum1(moles)
        mole_fractions = moles/tot_moles
        return mole_fractions

    def partial_pressures(self, w_i, p):
        # Methode to convert mass fractions into partial pressures
        return CasADi.SX(self.moleFractions(w_i)*p)

    def concentrations(self, w_i, T, p):
        # Methode to convert mass fractions into concentrations
        p_i = self.partial_pressures(w_i, p)
        return p_i/(self.R*T)

    def rho_cp_eff(self,w_i, T, p):
        # Methode to calculate the effective heat capacity of the cat + fluid phase
        # Eq. (32)
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)
        rho_fl = self.rho_fl(w_i, T, p)

        rho_cat = self.cat.get_property(Component.DENSITY, T)
        cp_cat = self.cat.get_property(Component.HEAT_CAPACITY, T)

        rho_cp_eff = (1-self.eps) * rho_cat * cp_cat + self.eps * rho_fl * cp_fl
        return rho_cp_eff

    def Re(self, w_i, T, u, p):
        # Methode for Reynolds number
        # Eq. (37)
        eta_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        d_p = self.cat_diameter
        rho_fl = self.rho_fl(w_i, T, p)

        Re = u * d_p * rho_fl / eta_fl
        return Re

    def Re_0 (self, T, p, u, w_i):
        # Methode for Reynolds number with velocity corrected for emtpy space velocity
        # Eq. (43)
        eps = self.eps
        eta_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        d_p = self.cat_diameter
        rho_fl = self.rho_fl(w_i, T, p)

        u_0 = u * eps

        Re_0 = u_0 * d_p * rho_fl / eta_fl
        return Re_0

    def Pe_0 (self, T, p, u, w_i):
        # Methode for Peclet number with velocity corrected for empty space velocity for thermal conductivity
        # Eq. (40)
        eps = self.eps
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)
        lambda_fl = self.calc_fluidThermalConductivity(T, w_i)
        d_p = self.cat_diameter
        rho_fl = self.rho_fl(w_i, T, p)

        u_0 = u * eps

        Pe_0 = u_0 * rho_fl * cp_fl * d_p / lambda_fl
        return Pe_0

    def Pe_0_diff(self, u, diffCoff):
        # Methode for Peclet number with velocity corrected for empty space velocity for diffusion
        # Eq. (17)
        eps = self.eps
        d_p = self.cat_diameter

        u_0 = u * eps

        Pe_0_diff = u_0 * d_p / diffCoff
        return Pe_0_diff

    def Pr(self, w_i, T):
        # Methode for Prandtl number
        # Eq. (38)
        eta_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)
        lambda_fl = self.calc_fluidThermalConductivity(T, w_i)

        Pr = eta_fl * cp_fl / lambda_fl
        return Pr

    def calc_fluidThermalConductivity(self, T, w_i):
        # Methode to calculate the thermal conductivity of the fluid phase with utilization viscosity mixing rule
        # Eq. (52)
        # Get component properties
        viscosity = []
        thermal_conductivity = []
        molar_weight = []
        for component in self.components:
            viscosity.append(component.get_property(Component.DYNAMIC_VISCOSITY, T))
            thermal_conductivity.append(component.get_property(Component.THERMAL_CONDUCTIVITY, T))
            molar_weight.append(component.get_property(Component.MOLECULAR_WEIGHT))
        moleFraction = self.moleFractions(w_i)

        # Calculate the gas mixture thermal conductivity employing viscosity mixing rule
        fluid_thermal_conductivity = (sum(moleFraction[comp_i] * thermal_conductivity[comp_i] / (
            sum(moleFraction[comp_j] *
                self.__calc_compActivity(molar_weight[comp_i], molar_weight[comp_j], viscosity[comp_i],
                                         viscosity[comp_j])
                for comp_j in range(len(self.components))))
                                  for comp_i in range(len(self.components))))
        return fluid_thermal_conductivity

    def __calc_compActivity(self, molar_weight_i, molar_weight_j, viscosity_i, viscosity_j):
        # Component activity for viscous mixing rule utilized in calculating thermal conductivity of a fluid mixture
        # Eq. (53)
        compActivity = (
                (1 + CasADi.sqrt(viscosity_i / viscosity_j) * (molar_weight_j / molar_weight_i) ** (1 / 4)) ** 2 /
                (CasADi.sqrt(8 * (1 + molar_weight_i / molar_weight_j))))
        return compActivity



