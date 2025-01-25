# Class to calculate molecular diffusion coefficient with Fuller equation
# and mixture averaging
# Till Kasselmann 14.01.2025
from V1.classes.ReactorSpecificQuantities.Component.Component import Component


class MixtureAveragedDiffusionCoefficient:

    def __init__(self, log, RSQ, GCF):
        self.log = log
        self.RSQ = RSQ
        self.GCF = GCF

    def calc(self, w_i, T, p, i): # i = 0,1,.... depending on the component index
        # get diffusion volumes, molar masses and densities for all components
        diff_volumes = []
        molar_weights = []
        densities = []
        components = self.RSQ.getComponents()
        for component in components:
            diff_volumes.append(component.get_diffusion_volume())
            molar_weights.append(component.get_molecular_weight())
            densities.append(component.get_density(T))

        # get fluid mixture values for the density and molar mass
        rho_fl = self.GCF.rho_fl(w_i, T, p)
        Mw_fl = self.GCF.massFraction_weighted_average(w_i, Component.MOLECULAR_WEIGHT)

        # calculating sum of rho_j/ Mw_j * binary diffusion coefficients_i,j
        summ = 0
        for j in enumerate(components):
            if j != i:
                fuller_ij = self.calc_Fuller(T, p, molar_weights[i], molar_weights[j[0]], diff_volumes[i], diff_volumes[j[0]])
                summ += densities[j[0]] / (fuller_ij * molar_weights[j[0]])

        avgDiffusionCoefficient = ((rho_fl/Mw_fl)-(densities[i]/molar_weights[i]))/summ
        return avgDiffusionCoefficient



    def calc_Fuller(self, T, p,  molar_mass_i, molar_mass_j, diff_volume_i, diff_volume_j):
        # Conversion from SI units to equation units
        p = p * 1e-5                           # Pa -> bar
        molar_mass_i = molar_mass_i * 1e3      # kg/mol -> g/mol
        molar_mass_j = molar_mass_j * 1e3      # kg/mol -> g/mol

        # Diffusion coefficient of i in j after Fuller in m^2/s
        DiffCoff = (1e-4 * (0.00143 * T ** 1.75 * ((1/molar_mass_i) + (1/molar_mass_j)) ** 0.5) /
                    (p * 2 ** 0.5 * (diff_volume_i ** (1 / 3) + diff_volume_j ** (1 / 3)) ** 2))

        return DiffCoff

    