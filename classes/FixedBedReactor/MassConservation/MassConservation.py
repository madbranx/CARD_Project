class MassConservation:
    def __init__(self, log, dimension, RSQ):
        self.log = log
        self.log.addEntry("initializing mass conservation", 2)
        self.dimension = dimension
        self.RSQ = RSQ


    def createCasADi(self, ae, T, w_i, u, p):
        self.log.addEntry("creating CasADi mass conservation equations (AE)", 3)
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ae, T, w_i, u)
        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D(ae, T, w_i, u)

    def __createCasADi_1D(self, ae, T, w_i, u):

        T_in = self.RSQ.getParameterValue('T_initial')
        u_in = self.RSQ.getParameterValue('u_initial')
        w_i_in = self.RSQ.getParameterValue('w_i_initial')

        for z in range(ae.size()[0]):
            ae[z] = 1  # TODO finish when f_rho etc is implemented

    def __createCasADi_2D(self, ae, temperature, w_i, u):
        pass
