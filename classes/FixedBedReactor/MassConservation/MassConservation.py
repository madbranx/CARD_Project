class MassConservation:
    def __init__(self, log, dimension, RSQ, GCF):
        self.log = log
        self.log.addEntry("initializing mass conservation", 2)
        self.dimension = dimension
        self.RSQ = RSQ
        self.GCF = GCF


    def createCasADi(self, ae, T, w_i, u, p):
        self.log.addEntry("creating CasADi mass conservation equations (AE)", 3)
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ae, T, w_i, u, p)
        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D(ae, T, w_i, u, p)

    def __createCasADi_1D(self, ae, T, w_i, u, p):

        T_in = self.RSQ.getParameterValue('T_initial')
        w_i_in = self.RSQ.getParameterValue('w_i_initial')
        p_in = self.RSQ.getParameterValue('p_initial')
        rho_fl_in = self.GCF.rho_fl(w_i_in, T_in, p_in)

        u_in = self.RSQ.getParameterValue('u_initial')

        for z in range(ae.size()[0]):
            rho_fl = self.GCF.rho_fl(w_i[z, :].T, T[z], p[z])
            ae[z] = u[z] - u_in * rho_fl/rho_fl_in


    def __createCasADi_2D(self, ae, temperature, w_i, u, p):
        pass
        # TODO
