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
        lambda_fl = self.calc_fluidThermalConductivity(T, w_i)
        lambda_cat = self.cat.get_property(Component.THERMAL_CONDUCTIVITY, T)

        Re = self.Re(w_i, T, u, p)
        Pr = self.Pr(w_i, T)

        Pe_eff = 1/((1/Pe_z_fl)+((lambda_cat/lambda_fl)/(Re/Pr)))
        return Pe_eff


    ## REACTION HEAT
    def reactionHeat(self, T, w_i, p):
        eff_factor = self.species_conservation.effFactor(w_i, T, p)
        reaction_rate = self.rate_equation(w_i, T, p)

        return -(1-self.eps) * eff_factor * reaction_rate * self.reactionEnthalpy


    ## RADIAL HEAT CONDUCTION
    def radial_heatConduction(self, radial_discretization, r, u_center, wTpu, wTpu_in=None, wTpu_out=None):
        T = wTpu[1]

        # Get Radial Discretization Values
        radial_centroids = radial_discretization.get_centroids()
        radial_faces = radial_discretization.get_faces()
        radial_faces_diff = radial_discretization.get_differences_faces()
        radial_centroids_diff = radial_discretization.get_differences_centroids()

        r_centroid = radial_centroids[r]
        r_face_in = radial_faces[r]
        r_face_out = radial_faces[r + 1]
        diff_faces = radial_faces_diff[r]

        # Calculate Thermal Conductivity of the current cell
        lambda_radial_center = self.calc_effective_radial_thermal_conductivity(wTpu, u_center, r_centroid)

        if wTpu_in is None: # Symmetry Boundary Condition
            q_r_in = 0
        else:
            T_in = wTpu_in[1]
            # Calculate Thermal Conductivity of the inlet cell
            lambda_radial_inlet = self.calc_effective_radial_thermal_conductivity(wTpu_in, u_center, radial_centroids[r - 1])

            # Calculate the mean of the lambda inlet and current
            diff_faces_inlet = radial_faces_diff[r-1]
            diff_faces_center = radial_faces_diff[r]
            lambda_radial_in = (lambda_radial_inlet*diff_faces_inlet + lambda_radial_center*diff_faces_center) / (
                    diff_faces_inlet + diff_faces_center)

            # calculating heat flux over inlet face
            diff_centroids_in = radial_centroids_diff[r - 1]
            q_r_in = -lambda_radial_in * (T_in - T) / diff_centroids_in

        if wTpu_out is None: # Wall Boundary Condition
            alpha_wall = self.calc_heatTransferCoefficient_wall(wTpu)
            T_inner_wall = self.calc_innerWallTemperature(wTpu)
            q_r_out = -alpha_wall * (T - T_inner_wall)
        else:
            T_out = wTpu_out[1]
            # Calculate Thermal Conductivity of the outlet cell
            lambda_radial_outlet = self.calc_effective_radial_thermal_conductivity(wTpu_out, u_center, radial_centroids[r + 1])

            # Calculate the mean of the lambda outlet and current
            diff_faces_center = radial_faces_diff[r]
            diff_faces_outlet = radial_faces_diff[r+1]
            lambda_radial_out = (lambda_radial_outlet * diff_faces_outlet + lambda_radial_center * diff_faces_center) / (
                        diff_faces_outlet + diff_faces_center)

            # calculating heat flux over outlet face
            diff_centroids_out = radial_centroids_diff[r]
            q_r_out = -lambda_radial_out * (T - T_out)  / diff_centroids_out

        # calculating total heat conduction
        radial_heat_conduction = 1 / r_centroid * (q_r_in * r_face_in - q_r_out * r_face_out) / diff_faces
        return radial_heat_conduction


    ## Methods for calculating radial heat conduction

    def calc_effective_radial_thermal_conductivity (self, wTpu, u_center, r_centroid):
        [w_i, T, p, u] = wTpu
        reactor_radius = self.reactorDiameter/2
        cat_diameter = self.cat_diameter

        thermalConductivity_fluid = self.calc_fluidThermalConductivity(T, w_i)
        thermalConductivity_bed = self.calc_thermalConductivity_bed(T, w_i)


        Reynold = self.Re_0(T, p, u, w_i)
        Peclet = self.Pe_0(T, p, u, w_i)

        K2 = 0.44 + 4 * CasADi.exp(-Reynold/70)

        f =  CasADi.if_else((reactor_radius - r_centroid) <= K2 * cat_diameter,         # Condition
                                ((reactor_radius - r_centroid) / K2 * cat_diameter) ** 2,     # if True
                                1                                                             # if False
                            )

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

    def calc_heatTransferCoefficient_wall(self, wTpu):
        [w_i, T, p, u] = wTpu
        reactorDiameter = self.reactorDiameter
        catDiameter = self.cat_diameter
        lambda_fl = self.calc_fluidThermalConductivity(T, w_i)
        lambda_bed = self.calc_thermalConductivity_bed(T, w_i)
        Reynold = self.Re_0(T, p, u, w_i)
        Prandtl = self.Pr(w_i, T)

        Nusselt = (1.3 + 5 * catDiameter/ reactorDiameter) * lambda_bed / lambda_fl + 0.19 * Reynold**0.75 * Prandtl**(1/3)

        alpha_wall = Nusselt * lambda_fl / catDiameter
        return alpha_wall

    def calc_resistanceWall(self):
        k_jacket = 1 / (1 / self.reactor_thermalConductivity * CasADi.log(((self.reactorDiameter/2) + self.reactor_wallThickness) / (self.reactorDiameter/2)))
        return k_jacket

    def calc_innerWallTemperature(self, wTpu):
        T = wTpu[1]
        k_jac = self.calc_resistanceWall()
        alpha = self.calc_heatTransferCoefficient_wall(wTpu)
        s_jac = self.reactor_wallThickness
        T_wall_out = self.T_wall
        T_wall_in = (alpha * T + k_jac/s_jac * T_wall_out) / (alpha + k_jac/s_jac)
        return T_wall_in