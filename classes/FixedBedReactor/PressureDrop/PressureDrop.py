class PressureDrop:
    def __init__(self, log, dimension, RSQ):
        self.log = log
        self.dimension = dimension
        self.RSQ = RSQ

    def createCasADi(self, ae, T, w_i, u, p):
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D()
        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D()

    def __createCasADi_1D(self):
        pass

    def __createCasADi_2D(self):
        pass