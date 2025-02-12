from classes.Properties.Kinetics import Kinetics
import casadi as CasADi

class MassConservation(Kinetics):
    def __init__(self):
        super().__init__()

    def massConservation(self, T, w_i, p):
            rho_fl = self.rho_fl(w_i, T, p)
            rho_fl_in = self.rho_fl(self.w_i_in, self.T_in, self.p_in)

            return rho_fl_in/rho_fl
