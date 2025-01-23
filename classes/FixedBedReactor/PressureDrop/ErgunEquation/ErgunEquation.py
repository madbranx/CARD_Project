import numpy as np
from classes.ReactorSpecificQuantities.Component.Component import Component

class ErgunEquation:
    def __init__(self, log, GCF, RSQ):
        self.log = log
        self.GCF = GCF
        self.RSQ = RSQ

    # Following Methods use location specific Arguments -> can be used for 1D and 2D!
    def calc(self, T, w_i, u, p):
        eps = self.RSQ.getParameterValue("bed_void_fraction")
        d_cat = self.RSQ.getParameterValue("cat_diameter")

        A = self.__A(w_i, T, eps, d_cat)
        B = self.__B(w_i, T, p, eps, d_cat)

        return A * u + B * np.pow(u, 2)

    def __A(self,w_i, T, eps, d_cat):
        viscosity_fl = self.GCF.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY,  T)

        A = (np.pow(1 - eps, 2) / np.pow(eps, 3)) * (150 * viscosity_fl / np.pow(d_cat, 2))
        return A

    def __B(self, w_i, T, p, eps, d_cat):
        rho_fl = self.GCF.rho_fl(w_i, T, p)

        B = ((1 - eps) / (np.pow(eps, 3))) * (1.75 / d_cat) * rho_fl
        return B
