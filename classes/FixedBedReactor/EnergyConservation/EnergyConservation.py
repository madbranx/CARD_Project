from classes.FixedBedReactor.EnergyConservation.convectiveHeatFlux.convectiveHeatFlux import ConvectiveHeatFlux
from classes.FixedBedReactor.EnergyConservation.heatConduction.EffAxialThermalConduction.EffAxialThermalConductivity import \
    EffAxialThermalConductivity
from classes.FixedBedReactor.EnergyConservation.reactionHeat.reactionHeat import ReactionHeat


class EnergyConservation:
    def __init__(self, log, dimension, RSQ, GCF, disc_axial, disc_radial = None):
        self.log = log
        self.log.addEntry("initializing energy conservation", 2)
        self.dimension = dimension
        self.RSQ = RSQ
        self.GCF = GCF
        self.disc_axial = disc_axial
        self.disc_radial = disc_radial

        # Initialize Calculation Classes
        self.convectiveHeatFlux = ConvectiveHeatFlux(log, self.GCF)
        self.effAxialHeatConductivity = EffAxialThermalConductivity(log, self.GCF, self.RSQ)
        self.effRadialHeatConductivity = None
        self.reactionHeat = ReactionHeat(log, self.GCF, self.RSQ)

    def createCasADi(self, ode, T, w_i, u, p):
        self.log.addEntry("creating CasADi energy conservation equations (ODE)", 3)
        from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor

        if self.dimension == FixedBedReactor.ONE_D:
            self.__createCasADi_1D(ode, T, w_i, u, p)

        elif self.dimension == FixedBedReactor.TWO_D:
            self.__createCasADi_2D()

    def __createCasADi_1D(self, ode, T, w_i, u, p):
        delta_axial = self.disc_axial.get_differences()
        T_in = self.RSQ.getParameterValue("T_in")

        for z in range(ode.size()[0]):
            rho_cp_eff = self.GCF.rho_cp_eff(w_i[z, :].T, T[z], p[z])

            if z == 0:  # Boundary Condition
                left_side = (T[z] - T_in) * rho_cp_eff
                conv_HF = (T[z] - T_in) /delta_axial[z] * self.convectiveHeatFlux.calc(T[z], w_i[z,:].T, u[z], p[z])
                axial_heatConduction = (T[z] - T_in) / pow(delta_axial[z],2) * (-self.effAxialHeatConductivity.calc(T[z], w_i[z,:].T, u[z], p[z]))
            else:
                left_side = (T[z] - T[z-1]) * rho_cp_eff
                conv_HF = (T[z] - T[z-1]) / delta_axial[z] * self.convectiveHeatFlux.calc(T[z], w_i[z,:].T, u[z], p[z])
                axial_heatConduction = (T[z] - T[z-1]) / pow(delta_axial[z], 2) * (-self.effAxialHeatConductivity.calc(T[z], w_i[z,:].T, u[z], p[z]))

            # 1D radial heat conduction using overall lambda
            # TODO correct?
            lambda_radial = self.RSQ.getParameterValue("lambda_total")
            d_reactor = self.RSQ.getParameterValue("reactorDiameter")
            T_wall = self.RSQ.getParameterValue("T_wall")
            radial_heatConduction = 4 * lambda_radial / d_reactor * (T[z]-T_wall)

            reactionHeat = self.reactionHeat.calc(T[z], w_i[z,:].T, p[z])


            ode[z] =    (
                        - left_side
                        - conv_HF
                        #- axial_heatConduction
                        #- radial_heatConduction
                        #- reactionHeat
                        )

    def __createCasADi_2D(self):
        pass