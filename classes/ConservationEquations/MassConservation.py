from classes.Properties.Kinetics import Kinetics
import casadi as CasADi

"""
The MassConservation class calculates the fluid density at the inlet and the current volume.
The factor of the fluid densitys is required for the AE of the mass conservation found in the FixedBedReactor class.
"""

class MassConservation(Kinetics):
    def __init__(self):
        super().__init__()

    def massConservation(self, T, w_i, p):
            rho_fl = self.rho_fl(w_i, T, p)
            rho_fl_in = self.rho_fl(self.w_i_in, self.T_in, self.p_in)

            return rho_fl_in/rho_fl
