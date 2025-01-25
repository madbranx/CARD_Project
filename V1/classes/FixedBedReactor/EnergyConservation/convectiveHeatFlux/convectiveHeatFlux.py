class ConvectiveHeatFlux:

    def __init__(self,log, GCF):
        self.log = log
        self.GCF = GCF

    # Following Methods use location specific Arguments -> can be used for 1D and 2D!
    def calc(self, T, w_i, u, p):
        from V1.classes.ReactorSpecificQuantities.Component.Component import Component
        rho_fl = self.GCF.rho_fl(w_i, T, p)
        cp_fl = self.GCF.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY,T)

        return u*rho_fl*cp_fl
