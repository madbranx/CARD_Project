# Class to calculate molecular diffusion coefficient with Fuller equation
# and mixture averaging
# Till Kasselmann 14.01.2025

class MolecularDiffCoff:

    def __init__(self):
        pass

    def calc_fuller_diffcoff(self, T, p,  molar_mass_i, molar_mass_j, diff_volume_i, diff_volume_j):
        # Conversion from SI units to equation units
        p = p * 1e-5                            # Pa -> bar
        molar_mass_i = molar_mass_i * 1e2      # kg/mol -> g/mol
        molar_mass_j = molar_mass_j * 1e2      # kg/mol -> g/mol

        # Diffusion coefficient of i in j after Fuller in m^2/s
        DiffCoff = ( 1e-4 * (0.00143 * T**1.75 * ((1/molar_mass_i) + (1/molar_mass_j))**(0.5)) /
                     (p * (2)**(0.5) * (diff_volume_i**(1/3) + diff_volume_j**(1/3))**2))

        return DiffCoff

    