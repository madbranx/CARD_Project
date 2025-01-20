class SpeciesConservation:
    def __init__(self, log, dimension, RSQ, GCV):
        self.log = log
        self.log.addEntry("initializing species conservation", 2)
        self.dimension = dimension
        self.RSQ = RSQ
        self.GCV = GCV

    def createCasADi(self, ode, T, w_i, u, p):
        self.log.addEntry("creating CasADi species conservation equations (ODE)", 3)
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ode)
        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D()

    def __createCasADi_1D(self, ode):
        for comp in range(4):
            for z in range(ode.size()[0]):
                ode[z, comp] = 1

    def __createCasADi_2D(self):
        pass