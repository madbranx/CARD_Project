from classes.ConservationEquations.SpeciesConservation import SpeciesConservation
from classes.Parameters.Component import Component
from classes.Properties.Kinetics import Kinetics


class EnergyConservation(Kinetics):
    def __init__(self):
        super().__init__()
        self.species_conservation = SpeciesConservation()


    ## CONVECTIVE HEAT FLUX
    def convectiveHeatFlux(self, T, w_i, u, p):
        rho_fl = self.rho_fl(w_i, T, p)
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)

        return u*rho_fl*cp_fl

    ## AXIAL HEAT CONDUCTION
    def effAxialThermalConductivity(self, T, w_i, u, p):
        rho_fl = self.rho_fl(w_i, T, p)
        cp_fl = self.massFraction_weighted_average(w_i, Component.HEAT_CAPACITY, T)

        Pe_eff = self.__Pe_eff(T, w_i, u, p)

        lambda_eff_axial = u + rho_fl * cp_fl * self.cat_diameter / Pe_eff
        return lambda_eff_axial

    def __Pe_eff(self, T, w_i, u, p):
        Pe_z_fl = 2
        lambda_fl = self.massFraction_weighted_average(w_i, Component.THERMAL_CONDUCTIVITY, T)
        lambda_cat = self.cat.get_thermal_conductivity(T)

        Re = self.Re(w_i, T, u, p)
        Pr = self.Pr(w_i, T)

        Pe_eff = 1/((1/Pe_z_fl)+((lambda_fl/lambda_cat)/(Re/Pr)))
        return Pe_eff

    ## RADIAL HEAT CONDUCTION
    def effRadialThermalConductivity(self, T, w_i, u, p):
        # TODO
        pass


    def wall_HeatTransferCoefficient(self, T, w_i, u, p):
        # TODO
        pass

    ## REACTION HEAT
    def reactionHeat(self, T, w_i, p):
        eff_factor = self.species_conservation.effFactor(w_i, T, p)
        reaction_rate = self.rate_equation(w_i, T, p)

        return (1-self.eps) * eff_factor * reaction_rate * self.reactionEnthalpy
