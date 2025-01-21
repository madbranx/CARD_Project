import numpy as np

from classes.ReactorSpecificQuantities.Component.Component import Component


class PressureDrop:
    def __init__(self, log, dimension, RSQ, GCF, disc_z, disc_r):
        self.log = log
        self.log.addEntry("initializing pressure drop", 2)
        self.dimension = dimension
        self.RSQ = RSQ
        self.GCF = GCF
        self.disc_z = disc_z
        self.disc_r = disc_r

    def createCasADi(self, ae, T, w_i, u, p):
        self.log.addEntry("creating CasADi pressure drop equations (AE)", 3)
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ae, T, w_i, u, p)
        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D(ae, T, w_i, u, p)

    def __createCasADi_1D(self, ae, T, w_i, u, p):
        self.__Ergun1D(ae, T, w_i, u, p)


    def __Ergun1D(self, ae, T, w_i, u, p):
        eps = self.RSQ.getParameterValue("bed_void_fraction")
        d_cat = self.RSQ.getParameterValue("cat_diameter")
        delta_z = self.disc_z.get_differences()
        for z in range(ae.size()[0]):

            A = self.__A(z, w_i, T, eps, d_cat)
            B = self.__B(z, w_i, T, p, eps, d_cat)

            factor = delta_z[z] * (A * u[z] + B * np.pow(u[z], 2))

            if z == 0: # Boundary condition
                p_in = self.RSQ.getParameterValue("p_in")
                ae[z] = (p[z] - p_in) + factor
            else:
                ae[z] = (p[z] - p[z-1]) + factor

    def __A(self,z, w_i, T, eps, d_cat):
        visc_fl = self.GCF.massFraction_weighted_average(w_i[z, :].T, Component.DYNAMIC_VISCOSITY,  T[z])

        A = (np.pow(1 - eps, 2) / np.pow(eps, 3)) * (150 * visc_fl / np.pow(d_cat, 2))
        return A

    def __B(self, z, w_i, T, p, eps, d_cat):
        rho_fl = self.GCF.rho_fl(w_i[z, :].T, T[z], p[z])

        B = ((1 - eps) / (np.pow(eps, 3))) * (1.75 / d_cat) * rho_fl
        return B

    def __createCasADi_2D(self, ae, T, w_i, u, p):
        pass  # TODO
