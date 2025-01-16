# Class to calculate molecular diffusion coefficient with Fuller equation
# and mixture averaging
# Till Kasselmann 14.01.2025

class MolecularDiffCoff:

    def __init__(self):
        pass

    def calc_fuller(self, T, p,  molar_mass_i, molar_mass_j, diff_volume_i, diff_volume_j):
        # T in K
        # p in Pa
        DiffCoff = (1.013e-7 * T**1.75 * ((molar_mass_i * 1e3)^-1 + (molar_mass_j * 1e3)^-1)**(-0.5) /
                    ( p * 10e-5 * (diff_volume_i**(1/3) + diff_volume_j**(1/3))))
        return DiffCoff

    