from classes.Parameters.Parameters import Parameters
from classes.Parameters.Component import Component

import casadi as CasADi


class FluidProperties(Parameters):

    def __init__(self):
        super().__init__()

    def rho_comp(self, w_i, T, p):
        return self.concentrations(w_i, T, p)*self.getMolarWeights()

    def rho_fl(self, w_i, T, p):
        Mw_fl = self.massFraction_weighted_average(w_i, Component.MOLECULAR_WEIGHT)
        rho = self.__c_fl(T, p) * Mw_fl
        return rho

    # mean fluid concentration
    def __c_fl(self, T, p):
        cf = p / (self.R * T)
        return cf

    def getMolarWeights(self):
        Mw_i = []
        for component in self.components:
            Mw_i.append(component.get_property(Component.MOLECULAR_WEIGHT))
        return Mw_i

    # mass fraction weighted average of a fluid property
    def massFraction_weighted_average(self, w_i, material_property, T=None):  # property = Component.DENSITY, ...
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
        Mw_i = self.getMolarWeights()
        moles = w_i / Mw_i
        tot_moles = CasADi.sum1(moles)
        mole_fractions = moles/tot_moles
        return mole_fractions

    def partial_pressures(self, w_i, p):
        return CasADi.SX(self.moleFractions(w_i)*p)

    def concentrations(self, w_i, T, p):
        p_i = self.partial_pressures(w_i, p)
        return p_i/(self.R*T)

    def rho_cp_eff(self,w_i, T, p):
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)
        rho_fl = self.rho_fl(w_i, T, p)

        rho_cat = self.cat.get_property(Component.DENSITY, T)
        cp_cat = self.cat.get_property(Component.HEAT_CAPACITY, T)

        rho_cp_eff = (1-self.eps) * rho_cat * cp_cat + self.eps * rho_fl * cp_fl
        return rho_cp_eff

    def Re(self, w_i, T, u, p):
        eta_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        d_p = self.cat_diameter
        rho_fl = self.rho_fl(w_i, T, p)

        Re = u * d_p * rho_fl / eta_fl
        return Re

    def Re_0 (self, T, p, u, w_i):
        # Reynolds number with velocity corrected for emtpy space velocity
        eps = self.eps
        eta_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        d_p = self.cat_diameter
        rho_fl = self.rho_fl(w_i, T, p)

        u_0 = u * eps

        Re_0 = u_0 * d_p * rho_fl / eta_fl
        return Re_0

    def Pe_0 (self, T, p, u, w_i):
        # Peclet number with velocity corrected for empty space velocity
        eps = self.eps
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)
        lambda_fl = self.calc_fluidThermalConductivity(T, w_i)
        d_p = self.cat_diameter
        rho_fl = self.rho_fl(w_i, T, p)

        u_0 = u * eps

        Pe_0 = u_0 * rho_fl * cp_fl * d_p / lambda_fl
        return Pe_0

    def Pr(self, w_i, T):
        eta_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)
        lambda_fl = self.calc_fluidThermalConductivity(T, w_i)

        Pr = eta_fl * cp_fl / lambda_fl
        return Pr

    def calc_fluidThermalConductivity(self, T, w_i):
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
        compActivity = (
                (1 + CasADi.sqrt(viscosity_i / viscosity_j) * (molar_weight_j / molar_weight_i) ** (1 / 4)) ** 2 /
                (CasADi.sqrt(8 * (1 + molar_weight_i / molar_weight_j))))
        return compActivity



