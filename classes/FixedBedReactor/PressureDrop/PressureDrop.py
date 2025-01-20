class PressureDrop:
    def __init__(self, log, dimension, RSQ, GCF):
        self.log = log
        self.log.addEntry("initializing pressure drop", 2)
        self.dimension = dimension
        self.RSQ = RSQ
        self.GCF = GCF

    def createCasADi(self, ae, T, w_i, u, p):
        self.log.addEntry("creating CasADi pressure drop equations (AE)", 3)
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ae)
        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D()

    def __createCasADi_1D(self, ae):
        for z in range(ae.size()[0]):
            ae[z] = 1

    def __createCasADi_2D(self):
        pass