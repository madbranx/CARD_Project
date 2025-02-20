from classes.Parameters.Component import Component
from classes.Properties.Kinetics import Kinetics

import casadi as CasADi

"""
The SpeciesConservation class calculates relevant parameters for the species conservation in the ode's.
The relevant equations and sources can be found in chapter 3 of the reactor equations overview.pdf file.
The effectiveness factor and thiele modulus are also calculated in this class. The kinetics are in their own Kinetics class.
"""

class SpeciesConservation(Kinetics):
    def __init__(self):
        super().__init__()

    ## AXIAL MASS FLOW
    def axialMassFlow(self, T, w_i, w_i_in, u, p, comp):
        # Calculate convective species flux
        rho_fl = self.rho_fl(w_i, T, p)

        # Eq. (11)
        j_i_ax = -u*rho_fl*(w_i_in[comp]-w_i[comp])
        return j_i_ax

    ## RADIAL MASS FLOW
    # Only Diffusion and Crossmixing
    # In the current state the radial dispersion leads to numerical problems when solving and is therefore currently commented out in FixedbedReactor.py

    def radialMassFlow(self, radial_discretization, r, comp, u_center, wTpu, wTpu_in=None, wTpu_out=None):
        # Calculate radial mass flow with two central differences (eq. (62-64))
        # Boundary conditions are considered by setting wTpu_(in/out) as None in the FixedBedReactor.py class.
        # If None value is given, the Neumann boundary condition are achieved by setting q_(in/out) = 0

        radial_centroids = radial_discretization.get_centroids()
        radial_faces = radial_discretization.get_faces()
        radial_faces_diff = radial_discretization.get_differences_faces()
        radial_centroids_diff = radial_discretization.get_differences_centroids()

        r_centroid = radial_centroids[r]
        r_face_in = radial_faces[r]
        r_face_out = radial_faces[r + 1]
        diff_faces = radial_faces_diff[r]


        if wTpu_in is None: # Symmetry Boundary Condition
            j_r_in = 0
        else:
            # Mass flow into the face of the center volume
            diff_centroids_in = radial_centroids_diff[r - 1]
            diff_faces_l = radial_faces_diff[r-1]
            diff_faces_r = radial_faces_diff[r]
            j_r_in = self.calc_j_dispersion(diff_centroids_in, diff_faces_l, diff_faces_r, wTpu_in, wTpu, comp, u_center, radial_centroids[r-1], r_centroid)

        if wTpu_out is None: # Wall Boundary Condition
            j_r_out = 0
        else:
            # Mass flow out of the face of the center volume
            diff_centroids_out = radial_centroids_diff[r]
            diff_faces_l = radial_faces_diff[r]
            diff_faces_r = radial_faces_diff[r+1]
            j_r_out = self.calc_j_dispersion(diff_centroids_out, diff_faces_l, diff_faces_r, wTpu, wTpu_out, comp, u_center, r_centroid, radial_centroids[r+1])

        # Calculate radial mass flow
        radial_mass_flow = 1 / r_centroid * (j_r_in * r_face_in - j_r_out * r_face_out) / diff_faces
        return radial_mass_flow

    def calc_j_dispersion(self, diff_centroids, diff_faces_l, diff_faces_r, wTpu_left, wTpu_right, comp, u_center, r_centroid_l, r_centroid_r):
        # Methode to calculate the mass flow between two volumes, refered to as left and right.
        [w_i_l, T_l, p_l, u_l] = wTpu_left
        [w_i_r, T_r, p_r, u_r] = wTpu_right

        # Calculate j_out
        rho_gas_l = self.rho_fl(w_i_l, T_l, p_l)
        rho_gas_r = self.rho_fl(w_i_r, T_r, p_r)
        # Calculate volume width weighted dispersion Coeff of left & right cell
        eff_DispCoff_l = self.calc_eff_disp_coff(comp, wTpu_left, u_center, r_centroid_l)
        eff_DispCoff_r = self.calc_eff_disp_coff(comp, wTpu_right, u_center, r_centroid_r)
        eff_DispCoff = (eff_DispCoff_l*diff_faces_l + eff_DispCoff_r*diff_faces_r)/(diff_faces_l+diff_faces_r)
        # Eq. (14)
        j = -eff_DispCoff * (rho_gas_l * w_i_l[comp] - rho_gas_r * w_i_r[comp]) / diff_centroids
        return j

    def calc_eff_disp_coff(self, comp, wTpu, u_center, r_centroid):
        # Methode to calculate the effective radial dispersion coefficient based on a model of Winterberg et al. .
        # Eq. (15-20)

        [w_i, T, p, u] = wTpu
        cat_diameter = self.cat_diameter
        reactor_radius = self.reactorDiameter / 2
        void_fraction = self.eps
        # Eq. (21+22)
        mix_DiffCoff = self.MixtureAveragedDiffusionCoefficient(w_i, T, p, comp)

        # Methode to calculate the radial thermal conductivity of the bed with stationary fluid based on a model of
        # Zehner, Bauer and Schlünder with eq. (20)
        bed_DiffCoff = (1 - CasADi.sqrt(1 - void_fraction)) * mix_DiffCoff

        # Eq. (17)
        Peclet = self.Pe_0_diff(u, mix_DiffCoff)
        Peclet_center = self.Pe_0_diff(u_center, mix_DiffCoff)

        # Eq. (16)
        K1 = 1 / 8 * 1/ (1 + 3/(Peclet_center**0.5))
        # Eq. (19)
        K2 = 0.44

        # Eq. (18)
        f = CasADi.if_else((reactor_radius - r_centroid) <= K2 * cat_diameter,  # Condition
                           ((reactor_radius - r_centroid) / (K2 * cat_diameter)) ** 2,  # if True
                           1  # if False
                           )

        # Eq. (15)
        rad_eff_disp_coeff = bed_DiffCoff + K1 * Peclet * u/u_center * f * mix_DiffCoff
        return rad_eff_disp_coeff

    def calc_sum_j(self, radial_discretization, r, u_center, wTpu, wTpu_in=None, wTpu_out=None):
        # Methode to calculate the sum of the radial mass flow of the volume for all species.
        # Sum is used to correct the mass conservation
        radial_faces_diff = radial_discretization.get_differences_faces()
        radial_centroids_diff = radial_discretization.get_differences_centroids()
        radial_centroids = radial_discretization.get_centroids()

        r_centroid = radial_centroids[r]

        sum_j = 0
        # Iterate through all components and sum the mass flow going in and out of the volume
        for comp in range(len(self.components)):
            if wTpu_in is None:  # Symmetry Boundary Condition
                j_r_in = 0
            else:
                diff_centroids_in = radial_centroids_diff[r - 1]
                diff_faces_l = radial_faces_diff[r - 1]
                diff_faces_r = radial_faces_diff[r]
                # From radial dispersion methode
                j_r_in = self.calc_j_dispersion(diff_centroids_in, diff_faces_l, diff_faces_r, wTpu_in, wTpu, comp, u_center, radial_centroids[r-1], r_centroid)

            if wTpu_out is None:  # Wall Boundary Condition
                j_r_out = 0
            else:
                diff_centroids_out = radial_centroids_diff[r]
                diff_faces_l = radial_faces_diff[r]
                diff_faces_r = radial_faces_diff[r + 1]
                # From radial dispersion methode
                j_r_out = self.calc_j_dispersion(diff_centroids_out, diff_faces_l, diff_faces_r, wTpu, wTpu_out, comp, u_center, r_centroid, radial_centroids[r+1])

            sum_j += (j_r_out - j_r_in)
        return sum_j

    ## CHANGE BY REACTION
    def changeByReaction(self, T, w_i, p, comp):
        # Shown in eq. (5)
        Mw_i = self.getMolarWeights()
        eff_factor = self.effFactor(w_i, T, p)
        reaction_rate = self.rate_equation(w_i, T, p)

        return (1 - self.eps) * Mw_i[comp] * self.nu[comp] * eff_factor * reaction_rate


    def effFactor(self, w_i, T, p):
        # Calculate effectiveness factor only for species CO2 as key component
        # Index of CO2 = 2
        comp = 2
        thiele = self.__calc_thiele(T, w_i, p, comp)
        # Eq. (27)
        effectiveness_factor = 3 / thiele * (1 / CasADi.tanh(thiele) - 1 / thiele)
        return effectiveness_factor

    def __calc_thiele(self, T, w_i, p, comp):
        # Methode to calculate thiele modulus only for species CO2 as key component
        # Thiele modulus with first order reaction
        c_i = self.concentrations(w_i, T, p)
        concentration_CO2 = c_i[2]
        stoichiometry_CO2 = self.nu[2]
        diameter_particle = self.cat_diameter
        reaction_rate = self.rate_equation(w_i, T, p)

        # Calculate the effective diffusion coefficient in porous cat with eq. (29)
        eff_diff_coff = self.calc_eff_diff_coff(T, w_i, p, comp)

        # Eq. (28)
        thiele = diameter_particle / 2 * ((stoichiometry_CO2 * reaction_rate) / (eff_diff_coff * concentration_CO2)) ** 0.5
        return thiele

    def calc_eff_diff_coff(self, T, w_i, p, comp):
        # Methode to calculate effective diffusion coefficient in porous cat with molecular and knudsen diffusion
        tortuosity_particle = self.cat_tortuosity
        porosity_particle = self.cat_porosity
        molar_mass = self.getMolarWeights()[comp]
        # With eq. (30)
        knudsen_diff_coff = self.Knudsen_diff_coff(T, molar_mass)
        # With eq. (21)
        molar_diff_coff = self.MixtureAveragedDiffusionCoefficient(w_i, T, p, comp)
        # Eq. (29)
        eff_diff_coff = 1 / ((tortuosity_particle ** 2 / porosity_particle) * (
                    1 / molar_diff_coff + 1 / knudsen_diff_coff))
        return eff_diff_coff

    def Knudsen_diff_coff(self, T, molarMass):
        # Methode to calculate knudsen diffusion coefficient
        # Eq. (30)
        knudsen_diff_coff = self.diameter_pore / 3 * ((8 * self.R * T) / (self.pi * molarMass)) ** 0.5
        return knudsen_diff_coff

    def MixtureAveragedDiffusionCoefficient(self, w_i, T, p, comp): # i = 0,1,.... depending on the component index
        # Methode to calculate mixture averaged diffusion coefficient based on the binary diffusion coefficient with
        # Fuller equation (eq. (22)). The mixture averaged diffusion coefficient is dependent on the component index.
        diff_volumes = []
        molar_weights = []
        densities = []

        for component in self.components:
            diff_volumes.append(component.get_property(Component.DIFFUSION_VOLUME))
            molar_weights.append(component.get_property(Component.MOLECULAR_WEIGHT))

        densities = self.rho_comp(w_i, T, p) # Using ideal Gas Law

        # Get fluid mixture values for the density and molar mass
        rho_fl = self.rho_fl(w_i, T, p)
        Mw_fl = self.massFraction_weighted_average(w_i, Component.MOLECULAR_WEIGHT)

        # Calculate sum of rho_j/ Mw_j * binary diffusion coefficients_i,j
        summ = 0
        for j in range(len(self.components)):
            if j != comp:
                fuller_ij = self.Fuller(T, p, molar_weights[comp], molar_weights[j], diff_volumes[comp], diff_volumes[j])
                summ += densities[j] / (fuller_ij * molar_weights[j])
        # Eq. (21)
        avgDiffusionCoefficient = ((rho_fl/Mw_fl) - (densities[comp] / molar_weights[comp])) / summ
        return avgDiffusionCoefficient

    def Fuller(self, T, p, molar_mass_i, molar_mass_j, diff_volume_i, diff_volume_j):
        # Methode to calculate binary diffusion coefficient with fuller equation.
        # Conversion from SI units to equation units
        p = p * 1e-5                           # Pa -> bar
        molar_mass_i = molar_mass_i * 1e3      # kg/mol -> g/mol
        molar_mass_j = molar_mass_j * 1e3      # kg/mol -> g/mol

        # Diffusion coefficient of i in j after Fuller in m^2/s
        # Eq. (22)
        DiffCoff = (1e-4 * (0.00143 * T ** 1.75 * ((1/molar_mass_i) + (1/molar_mass_j)) ** 0.5) /
                    (p * 2 ** 0.5 * (diff_volume_i ** (1 / 3) + diff_volume_j ** (1 / 3)) ** 2))
        return DiffCoff
