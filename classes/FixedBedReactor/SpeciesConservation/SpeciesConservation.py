from classes.FixedBedReactor.SpeciesConservation.AxialMassFlow.AxialMassFlow import AxialMassFlow
from classes.FixedBedReactor.SpeciesConservation.ChangeByReaction.ChangeByReaction import ChangeByReaction
from classes.FixedBedReactor.SpeciesConservation.RadialMassFlow.RadialMassFlow import RadialMassFlow


class SpeciesConservation:
    def __init__(self, log, dimension, RSQ, GCF, disc_axial, disc_radial = None):
        self.log = log
        self.log.addEntry("initializing species conservation", 2)
        self.dimension = dimension
        self.RSQ = RSQ
        self.GCF = GCF
        self.disc_axial = disc_axial
        self.disc_radial = disc_radial

        # Initialize Calculation Classes
        self.axialMassFlow = AxialMassFlow(log, self.GCF)
        self.radialMassFlow = RadialMassFlow(log, self.RSQ, self.GCF)
        self.changeByReaction = ChangeByReaction(log, self.RSQ, self.GCF)

    def createCasADi(self, ode, T, w_i, u, p):
        self.log.addEntry("creating CasADi species conservation equations (ODE)", 3)
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor

        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ode, T, w_i, u, p)

        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D()

    def __createCasADi_1D(self, ode, T, w_i, u, p):
        eps = self.RSQ.getParameterValue("bed_void_fraction")
        w_i_in = self.RSQ.getParameterValue("w_i_in")
        delta_axial = self.disc_axial.get_differences()

        n_components = self.RSQ.getNComponents()
        for comp in range(n_components):
            for z in range(ode.size()[0]):
                rho_fl = self.GCF.rho_fl(w_i[z, :].T, T[z], p[z])

                if z == 0:  # Boundary Condition
                    axialMassFlow = self.axialMassFlow.calc(T[z], w_i[z,:].T, w_i_in, u[z], p[z], comp)
                else:
                    axialMassFlow = self.axialMassFlow.calc(T[z], w_i[z, :].T, w_i[z-1,:].T, u[z], p[z], comp)

                ode[z, comp] = (
                                - axialMassFlow/ (delta_axial[z] * eps * rho_fl)
                                + self.changeByReaction.calc(T[z], w_i[z,:].T, p[z], comp) / (eps * rho_fl)
                                )

    def __createCasADi_2D(self):
        # TODO
        self.radialMassFlow.calc()
