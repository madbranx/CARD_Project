from classes.ConservationEquations.SpeciesConservation import SpeciesConservation
from classes.Parameters.Component import Component
from classes.Properties.Kinetics import Kinetics

import casadi as CasADi


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

        lambda_eff_axial = u * rho_fl * cp_fl * self.cat_diameter / Pe_eff
        return lambda_eff_axial

    def __Pe_eff(self, T, w_i, u, p):
        Pe_z_fl = 2
        lambda_fl = self.massFraction_weighted_average(w_i, Component.THERMAL_CONDUCTIVITY, T)
        lambda_cat = self.cat.get_property(Component.THERMAL_CONDUCTIVITY, T)

        Re = self.Re(w_i, T, u, p)
        Pr = self.Pr(w_i, T)

        Pe_eff = 1/((1/Pe_z_fl)+((lambda_cat/lambda_fl)/(Re/Pr)))
        return Pe_eff

    ## RADIAL HEAT CONDUCTION
    def effRadialThermalConductivity(self, T, T_in, T_out, p, u, u_center, w_i, delta_r_centeroids_in, delta_r_centeroids_out, delta_r_faces, r_face_in, r_face_out, r_centeroid):
        thermalConductivity_rad = self.calc_radialEffectiveThermalConductivity(T, p, u, u_center, w_i, r_centeroid)

        q_r_in = -thermalConductivity_radial * (T - T_in) / delta_r_centeroids_in
        q_r_out = -thermalConductivity_radial * (T_out - T) / delta_r_centeroids_out

        radial_heat_conduction = 1 / r_centeroid * (q_r_out * r_face_out - q_r_in * r_face_in) / delta_r_faces
        return radial_heat_conduction


    def wallRadialThermalConductivity(self, T, T_in, p, u, w_i, delta_r_centroids_in, delta_r_faces, r_face_in, r_face_out, r_centeroid):
        thermalConductivity_radial = self.calc_thermalConductivity_bed(T, w_i, p)
        alpha_wall = self.calc_heatTransferCoefficient_wall(T, p, u, w_i)
        T_inner_wall = self.calc_innerWallTemperature(T, p, u, w_i)

        q_r_in = -thermalConductivity_radial * (T - T_in) / delta_r_centroids_in
        q_r_out = -alpha_wall * (T_inner_wall - T)                   # TODO stimmt das so?

        radial_heat_conduction = 1 / r_centeroid * (q_r_out * r_face_out - q_r_in * r_face_in) / delta_r_faces
        return radial_heat_conduction




    ## REACTION HEAT
    def reactionHeat(self, T, w_i, p):
        eff_factor = self.species_conservation.effFactor(w_i, T, p)
        reaction_rate = self.rate_equation(w_i, T, p)

        return -(1-self.eps) * eff_factor * reaction_rate * self.reactionEnthalpy


    ## Methodes for calculating radial heat conduction

    def calc_radialEffectiveThermalConductivity (self, T, p, u, u_center, w_i, r_centeroid):
        reactor_radius = self.reactorDiameter/2
        cat_diameter = self.cat_diameter

        thermalConductivity_fluid = self.calc_fluidThermalConductivity(T, w_i)
        thermalConductivity_bed = self.calc_thermalConductivity_bed(T, w_i)


        Reynold = self.Re_0(T, p, u, w_i)
        Peclet = self.Pe_0(T, p, u, w_i)

        K2 = 0.44 + 4 * CasADi.exp(-Reynold/70)


        if 0 < (reactor_radius - r_centeroid) <= K2 * cat_diameter:
            f = ((reactor_radius - r_centeroid) / K2 * cat_diameter) ** 2
        elif K2 * cat_diameter < reactor_radius - r_centeroid <= reactor_radius:
            f = 1
        else:
            f = 1               # TODO error hinzufügen

        rad_eff_thermal_conductivity = thermalConductivity_bed + 1/8 * Peclet * u/u_center * f * thermalConductivity_fluid
        return rad_eff_thermal_conductivity

    def calc_thermalConductivity_bed (self, T, w_i):
        cat_diameter = self.cat_diameter
        eps = self.eps

        thermalConductivity_fluid = self.calc_fluidThermalConductivity(T, w_i)
        thermalConductivity_cat = self.cat.get_property(Component.THERMAL_CONDUCTIVITY, T)

        k_rad = 4 * self.radiation_blackBody / (2/self.cat_emissionCoefficient -1) * T**3 * cat_diameter / thermalConductivity_fluid
        k_p = thermalConductivity_cat / thermalConductivity_fluid
        B = 1.25 * ((1-eps) / eps) ** (10/9)
        N = 1 - B/k_p
        k_c = 2/N * (B/(N**2) * (k_p - 1)/k_p * CasADi.log(k_p/B) - (B + 1)/2 - (B - 1)/N )

        k_bed_stationary = 1 - CasADi.sqrt(1 - eps) + CasADi.sqrt(1 - eps) * k_c

        thermalConductivity_bed = (k_bed_stationary + (1 - CasADi.sqrt(1 - eps)) * k_rad + CasADi.sqrt(1 - eps) * 1 / (1/k_rad + 1/k_p)) * thermalConductivity_fluid
        return thermalConductivity_bed


    ## Methodes for calculating wall heat transfer coefficient with Nusselt correlation

    def calc_heatTransferCoefficient_wall(self, T, p, u, w_i):
        reactorDiameter = self.reactorDiameter
        catDiameter = self.cat_diameter
        lambda_fl = self.calc_fluidThermalConductivity(T, w_i)
        lambda_bed = self.calc_thermalConductivity_bed(T, w_i, p)
        Reynold = self.Re_0(T, p, u, w_i)
        Prandtl = self.Pr(w_i, T)

        Nusselt = (1.3 + 5 * reactorDiameter / catDiameter) * lambda_bed / lambda_fl + 0.19 * Reynold**0.75 * Prandtl**(1/3)

        alpha_wall = Nusselt * lambda_fl / lambda_bed
        return alpha_wall

    def calc_resistanceWall(self):
        k_jacket = 1 / (1 / self.reactor_thermalConductivity * CasADi.log(((self.reactorDiameter/2) + self.reactor_wallThickness) / (self.reactorDiameter/2)))
        return k_jacket

    def calc_innerWallTemperature(self, T, p, u, w_i):
        alpha_bed = self.calc_heatTransferCoefficient_wall(T, p, u, w_i)
        alpha_wall = self.calc_resistanceWall() / self.reactor_wallThickness
        T_wall = self.T_wall
        #return (alpha_bed * T + alpha_wall * T_wall) / (alpha_bed + alpha_wall)    #TODO stimmt nicht!
        return self.T_wall