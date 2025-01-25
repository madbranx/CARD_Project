from classes.Parameters.Component import Component
from classes.Properties.Kinetics import Kinetics

class PressureDrop(Kinetics):
    def __init__(self):
        super().__init__()

    def pressureDrop(self, T, w_i, u, p):
        # according to ERGUN
        A = self.__A(w_i, T)
        B = self.__B(w_i, T, p)

        return A * u + B * u**2

    def __A(self, w_i, T):
        viscosity_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        A = ((1 - self.eps)** 2) / self.eps**3 * (150 * viscosity_fl / self.cat_diameter ** 2)
        return A

    def __B(self, w_i, T, p):
        rho_fl = self.rho_fl(w_i, T, p)
        B = ((1 - self.eps) / self.eps**3) * (1.75 / self.cat_diameter) * rho_fl
        return B


