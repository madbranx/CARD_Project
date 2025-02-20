from classes.Parameters.Component import Component
from classes.Properties.Kinetics import Kinetics
import casadi as CasADi

"""
The PressureDrop class calculates the pressure drop with the ergun equation.
"""

class PressureDrop(Kinetics):
    def __init__(self):
        super().__init__()

    def pressureDrop(self, T, w_i, u, p):
        # Ergun equation (eq. (60)) split into viscous term (A) and turbulent term (B)
        A = self.__A(w_i, T)
        B = self.__B(w_i, T, p)

        return A * u + B * u**2

    def __A(self, w_i, T):
        # Viscous term of ergun equation
        viscosity_fl = self.massFraction_weighted_average(w_i, Component.DYNAMIC_VISCOSITY, T)
        A = ((1 - self.eps)** 2) / self.eps**3 * (150 * viscosity_fl / self.cat_diameter ** 2)
        return A

    def __B(self, w_i, T, p):
        # Turbulent term of ergun equation
        rho_fl = self.rho_fl(w_i, T, p)
        B = ((1 - self.eps) / self.eps**3) * (1.75 / self.cat_diameter) * rho_fl
        return B


