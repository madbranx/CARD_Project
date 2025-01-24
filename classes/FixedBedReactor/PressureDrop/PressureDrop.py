from classes.FixedBedReactor.PressureDrop.ErgunEquation.ErgunEquation import ErgunEquation


class PressureDrop:
    def __init__(self, log, dimension, RSQ, GCF, disc_axial, disc_radial):
        self.log = log
        self.log.addEntry("initializing pressure drop", 2)
        self.dimension = dimension
        self.RSQ = RSQ
        self.GCF = GCF
        self.disc_z = disc_axial
        self.disc_r = disc_radial

        # Initialize Calculation Classes
        self.PressureDrop = ErgunEquation(self.log, self.GCF, self.RSQ)

    def createCasADi(self, ae, T, w_i, u, p):
        self.log.addEntry("creating CasADi pressure drop equations (AE)", 3)
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ae, T, w_i, u, p)
        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D(ae, T, w_i, u, p)

    def __createCasADi_1D(self, ae, T, w_i, u, p):
        delta_z = self.disc_z.get_differences()

        for z in range(ae.size()[0]):
            pressureDrop = self.PressureDrop.calc(T[z], w_i[z, :].T, u[z], p[z])

            if z == 0:  # Boundary condition
                p_in = self.RSQ.getParameterValue("p_in")
                ae[z] = (p[z] - p_in) / delta_z[z] + pressureDrop
            else:
                ae[z] = (p[z] - p[z - 1]) / delta_z[z] + pressureDrop


    def __createCasADi_2D(self, ae, T, w_i, u, p):
        pass  # TODO
