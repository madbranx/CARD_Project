# Class to calculate axial mass flow
# Till Kasselmann 14.01.2025

class AxialMassFlow:

    def __init__(self,log, GCF):
        self.log = log
        self.GCF = GCF

    # Following Methods use location specific Arguments -> can be used for 1D and 2D!
    def calc(self, T, w_i, w_i_in, u, p, comp):
        # Only convection
        rho_fl = self.GCF.rho_fl(w_i, T, p)

        j_i_ax = -u*rho_fl*(w_i[comp]-w_i_in[comp])
        return j_i_ax

