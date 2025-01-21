import casadi as CasADi

from classes.ReactorSpecificQuantities.Component.Component import Component


class GeneralConversionFunctions:
    def __init__(self, log, RSQ):
        self.log = log
        self.RSQ = RSQ

    # mean fluid density
    def rho_fl(self, w_i, T, p):
        Mw_fl = self.massFraction_weighted_average(w_i, Component.MOLECULAR_WEIGHT)
        rho = self.__c_fl(T, p) * Mw_fl
        return rho

    # mean fluid concentration
    def __c_fl(self, T, p):
        R = self.RSQ.getParameterValue("R")
        cf = p / (R * T)
        return cf

    def massFraction_weighted_average(self, w_i, material_property, T=None):  # property = Component.DENSITY, ...
        components = self.RSQ.getComponents()
        property_values = []
        for component in components:
            property_values.append(component.get_property(material_property, T))

        unity = CasADi.DM.ones(len(property_values))
        return 1 / (CasADi.dot(w_i / CasADi.SX(property_values), unity))

