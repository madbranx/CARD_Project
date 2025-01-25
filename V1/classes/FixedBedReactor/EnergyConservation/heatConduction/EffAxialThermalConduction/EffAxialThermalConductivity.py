from V1.classes.ReactorSpecificQuantities.Component.Component import Component


class EffAxialThermalConductivity:

    def __init__(self,log, GCF, RSQ):
        self.log = log
        self.GCF = GCF
        self.RSQ = RSQ

    # Following Methods use location specific Arguments -> can be used for 1D and 2D!
    def calc(self, T, w_i, u, p):
        d_p = self.RSQ.getParameterValue("cat_diameter")
        rho_fl = self.GCF.rho_fl(w_i, T, p)
        cp_fl = self.GCF.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)

        Pe_eff = self.__Pe_eff(T, w_i, u, p)

        lambda_eff = u + rho_fl * cp_fl * d_p / Pe_eff
        return lambda_eff

    def __Pe_eff(self, T, w_i, u, p):
        Pe_z_fl = 2
        lambda_fl = self.GCF.massFraction_weighted_average(w_i, Component.THERMAL_CONDUCTIVITY, T)
        cat = self.RSQ.getCatalyst()
        lambda_cat = cat.get_thermal_conductivity(T)

        Re = self.GCF.Re(w_i, T, u, p)
        Pr = self.GCF.Pr(w_i, T)

        Pe_eff = 1/((1/Pe_z_fl)+((lambda_fl/lambda_cat)/(Re/Pr)))
        return Pe_eff