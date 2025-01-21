import casadi as CasADi

class GeneralConversionFunctions:
    def __init__(self, log, RSQ):
        self.log = log
        self.RSQ = RSQ

    # mean fluid density
    def rho_fl(self, w_i, T, p):
        rho = self.__c_fl(T, p) * self.__Mw_fl(w_i)
        return rho

    # mean fluid concentration
    def __c_fl(self, T, p):
        R = self.RSQ.getParameterValue("R")
        cf = p / (R * T)
        return cf

    # mean fluid molecular weight
    def __Mw_fl(self, w_i):
        Mw_i = []
        for component in self.RSQ.getComponents():
            Mw_i.append(component.get_molecular_weight())

        unity = CasADi.DM.ones(len(Mw_i))
        Mw_fl =  1 / (CasADi.dot(w_i / CasADi.SX(Mw_i), unity))
        return Mw_fl

    def massFraction_weighted_average(self, w_i, material_property, T=None):  # property = Component.DENSITY, ...
        components = self.RSQ.getComponents()

        summ = 0
        for component, w in zip(components, w_i):
            prop = component.get_property(material_property, T)
            summ = summ + (w / prop)

        return 1/summ

