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
    def effRadialThermalConductivity(self, T, T_in, T_out, p, w_i, delta_r_centeroids_in, delta_r_centeroids_out, delta_r_faces, r_face_in, r_face_out, r_centeroid):
        thermalConductivity_radial = self.calc_thermalConductivity_bed(T, w_i, p)

        q_r_in = -thermalConductivity_radial * (T - T_in) / delta_r_centeroids_in
        q_r_out = -thermalConductivity_radial * (T_out - T) / delta_r_centeroids_out

        radial_heat_conduction = 1 / r_centeroid * (q_r_out * r_face_out - q_r_in * r_face_in) / delta_r_faces
        return radial_heat_conduction


    def wallRadialThermalConductivity(self, T, T_in, p, w_i, delta_r_centroids_in, delta_r_faces, r_face_in, r_face_out, r_centeroid):
        thermalConductivity_radial = self.calc_thermalConductivity_bed(T, w_i, p)
        alpha_wall = self.calc_heatTransferCoefficient_contact(T, p, w_i)
        T_inner_wall = self.calc_innerWallTemperature(T, p, w_i)

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

    def calc_thermalConductivity_bed(self, T, w_i, p):
        k_G = self.calc_k_g(T, p, w_i)
        k_cat = self.calc_k_cat(T, w_i)
        k_rad = self.calc_k_rad(T, w_i)
        k_c = self.calc_k_c(T, p, w_i)

        thermalConductivity_bed = (((1 - CasADi.sqrt(1 - self.eps)) * self.eps * (1/(self.eps - 1 + 1/k_G) + k_rad)
                                    + CasADi.sqrt(1 - self.eps) * (self.cat_flatteningCoefficient * k_cat + (1 - self.cat_flatteningCoefficient) * k_c)
                                    ) * self.calc_fluidThermalConductivity(T, w_i))
        return thermalConductivity_bed

    def calc_k_c(self, T, p, w_i):
        N = self.__calc_N(T, p, w_i)
        B = self.__calc_deformationParameter()
        k_cat = self.calc_k_cat(T, w_i)
        k_rad = self.calc_k_rad(T, w_i)
        k_G = self.calc_k_g(T, p, w_i)

        term1 = ((2/N) * (
                (B * (k_cat + k_rad - 1)) / (N**2 * k_G * k_cat)
                * CasADi.log((k_cat + k_rad) / (B * (k_G + (1 - k_G) * (k_cat + k_rad))))
                ))
        term2 = (((B + 1) / (2 * B)) * (
                (k_rad / k_G) - B * (1 + ((1 - k_G) / k_G) * k_rad)
                ))
        term3 = (B - 1) / (N * k_G)
        return term1 + term2 - term3

    def __calc_N(self, T, p, w_i):
        N = (1 / self.calc_k_g(T, p, w_i) *
             (1 + (self.calc_k_rad(T, w_i) - self.__calc_deformationParameter() * self.calc_k_g(T, p, w_i)) / self.calc_k_cat(T, w_i))
             - self.__calc_deformationParameter() * (1 / self.calc_k_g(T, p, w_i) - 1) * (1 + self.calc_k_rad(T, w_i) / self.calc_k_cat(T, w_i)))
        return N

    def calc_k_g(self, T, w_i):
        return 1 / self.calc_fluidThermalConductivity(T, w_i)

    def __calc_meanFreePath(self, T, p, w_i):
        # Gas kinetics to calculate mean free path of gas mixture
        # effective collision area weighted by molar fraction
        collisionAreas = []
        moleFractions = self.moleFractions(w_i)

        for component in self.components:
            collisionAreas.append(component.get_property(Component.COLLISION_AREA))
        # mole fraction weighted effective collision area
        eff_collisionArea = sum(moleFractions[component] * collisionAreas[component] for component in range(len(self.components)))

        meanFreePath = 1 / CasADi.sqrt(2) * self.boltzmann * T / (p * eff_collisionArea)
        return meanFreePath

    def calc_k_cat(self, T, w_i):
        return self.cat.get_property(Component.THERMAL_CONDUCTIVITY) / self.calc_fluidThermalConductivity(T, w_i)

    def __calc_deformationParameter(self):
        return self.cat_shapeFactor * ((1 - self.eps)/self.eps) ** (10/9) * self.cat_distributionFunction

    def calc_k_rad(self, T, w_i):
        # Heat conduction parameter for radiation
        radiationCoefficient = self.radiation_blackBody
        emissionCoefficient = self.cat_emissionCoefficient
        particleDiameter = self.cat_diameter

        k_rad = 4 * radiationCoefficient / (2/emissionCoefficient - 1) * T**3 * particleDiameter / self.calc_fluidThermalConductivity(T, w_i)
        return k_rad

    ### Methodes for calculating wall heat transfer coefficient with Nusselt correlation

    def calc_heatTransferCoefficient_wall(self, T, p, u, w_i):
        reactorDiameter = self.reactorDiameter
        catDiameter = self.cat_diameter
        lambda_fl = self.calc_fluidThermalConductivity(T, w_i)
        lambda_bed = self.calc_thermalConductivity_bed(T, w_i, p)
        Reynold = self.Re(w_i, T, u, p)
        Prandtl = self.Pr(w_i, T)

        Nusselt = (1.3 + 5 * reactorDiameter / catDiameter) * lambda_bed/lambda_fl + 0.19 * Reynold**0.75 * Prandtl**(1/3)

        alpha_wall = Nusselt * lambda_fl / lambda_bed
        return alpha_wall

    def calc_resistanceWall(self):
        k_jacket = 1 / (1 / self.reactor_thermalConductivity * CasADi.log(((self.reactorDiameter/2) + self.reactor_wallThickness) / (self.reactorDiameter/2)))
        return k_jacket

    def calc_innerWallTemperature(self, T, p, w_i):
        alpha_bed = self.calc_heatTransferCoefficient_contact(T, p, w_i)
        alpha_wall = self.calc_resistanceWall() / self.reactor_wallThickness
        T_wall = self.T_wall
        #return (alpha_bed * T + alpha_wall * T_wall) / (alpha_bed + alpha_wall)    #TODO stimmt nicht!
        return self.T_wall

