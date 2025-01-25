from classes.Parameters.Parameters import Parameters
from classes.Parameters.Component import Component


class FluidProperties(Parameters):

    def __init__(self):
        super().__init__()


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
            summ += w_i[i] / property_value
        return 1/summ

    def moleFractions(self, w_i):
        Mw_i = self.getMolarWeights()
        Mw_fl = self.massFraction_weighted_average(w_i, Component.MOLECULAR_WEIGHT)
        return w_i*Mw_i/Mw_fl

    def partial_pressures(self, w_i, p):
        return self.moleFractions(w_i)*p

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

    def Pr(self, w_i, T):
        eta_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)
        lambda_fl = self.massFraction_weighted_average(w_i, Component.THERMAL_CONDUCTIVITY, T)

        Pr = eta_fl * cp_fl / lambda_fl
        return Pr
