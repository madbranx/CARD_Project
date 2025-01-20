from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor


class MassConservation:
    def __init__(self, log, dimension, RSQ):
        self.log = log
        self.dimension = dimension
        self.RSQ = RSQ

    def createCasADi(self, ae, temperature, w_i, u):
        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ae, temperature, w_i, u)
        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D(ae, temperature, w_i, u)

    def __createCasADi_1D(self, ae, T, w_i, u):

        T_in = self.RSQ.getParameter('u_in')

        for z in T.size()[0]:
            ae[z] = u[z] - u_in * f_rho(w_i_in, T_in) / f_rho(w_i[z, :].T, T[z])

    def __createCasADi_2D(self, ae, temperature, w_i, u):
        pass
