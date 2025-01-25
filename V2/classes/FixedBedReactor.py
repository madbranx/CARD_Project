import casadi as CasADi

from V2.classes.Parameters.Discretization import Discretization
from V2.classes.ConservationEquations.EnergyConservation import EnergyConservation
from V2.classes.ConservationEquations.MassConservation import MassConservation
from V2.classes.ConservationEquations.PressureDrop import PressureDrop
from V2.classes.ConservationEquations.SpeciesConservation import SpeciesConservation


class FixedBedReactor(EnergyConservation, MassConservation, PressureDrop, SpeciesConservation):
    def __init__(self, dimension, n_axial, n_radial=1):
        super().__init__()
        self.dimension = dimension

        # Discretizations
        if dimension == 1:
            n_radial = 1
        self.axial_discretization = Discretization(n_axial, start=0, end=self.reactorLength)
        self.radial_discretization = Discretization(n_radial, start=0, end=self.reactorDiameter/2)
        if self.dimension == 1:
            self.n_spatial = self.axial_discretization.num_volumes
        else:
            self.n_spatial = self.axial_discretization.num_volumes * self.radial_discretization.num_volumes
        self.n_components = len(self.components)

        # CasADi Variables
        self.w_i = None
        self.T = None
        self.p = None
        self.u = None

        # CasADi AEs/ODEs
        self.AE_m = None
        self.AE_p = None
        self.ODE_T = None
        self.ODE_wi = None

        # CasADi DAE
        self.DAE = None

        # SETUP
        self.__initializeCasADiStructure()
        self.__fillCasADiStructure()
        self.__reshapeCasADi_createODE()


    def __initializeCasADiStructure(self):
        # create CasADi Variables

        # Define differential variables
        self.w_i = CasADi.SX.sym('w_i', (self.n_spatial, self.n_components))
        self.T = CasADi.SX.sym('T', self.n_spatial)

        # Define algebraic variables
        self.u = CasADi.SX.sym('u', self.n_spatial)
        self.p = CasADi.SX.sym('p', self.n_spatial)

        # create CasADi AEs/ODEs
        self.ODE_wi = CasADi.SX.sym("ode_wi", (self.n_spatial, self.n_components))
        self.ODE_T = CasADi.SX.sym("ode_wi", self.n_spatial)
        self.AE_m = CasADi.SX.sym("ae_mass", self.n_spatial)
        self.AE_p = CasADi.SX.sym("ae_pressure_drop", self.n_spatial)

    def __fillCasADiStructure(self):
        w_i = self.w_i
        T = self.T
        p = self.p
        u = self.u

        radial_deltas = self.radial_discretization.get_differences()
        axial_deltas = self.axial_discretization.get_differences()

        for r, delta_r in enumerate(radial_deltas):
            for z, delta_z in enumerate(axial_deltas):
                current = z + r*len(axial_deltas)
                before_z = current -1
                before_r = current-len(axial_deltas)

            ## 1) MASS CONSERVATION
                self.AE_m[current] = u[current] - self.u_in * self.massConservation(T[current], w_i[current, :].T, p[current])

                ## 2) PRESSURE DROP
                if z==0:  # Inlet Boundary Condition
                    delta_p =  p[current] - self.p_in
                else:
                    delta_p =  p[current] - p[before_z]
                self.AE_p[current] = delta_p / delta_z + self.pressureDrop(T[current], w_i[current, :].T, u[current], p[current])

            ## 3) ENERGY CONSERVATION
                # 3.1) Axial Energy Conversation
                if z == 0:  # Inlet Boundary Condition
                    delta_T =  (T[current] - self.T_in)
                else:
                    delta_T =  (T[current] - T[before_z])

                axial_convectiveHeatFlux = self.convectiveHeatFlux(T[current], w_i[current, :].T, u[current], p[current]) * delta_T / delta_z

                lambda_eff_axial = self.effAxialThermalConductivity(T[current], w_i[current, :].T, u[current], p[current])
                axial_heatConduction = (-lambda_eff_axial) * delta_T / (delta_z ** 2)

                ## 3.2) Radial Energy Conversation
                if self.dimension == 2:
                    if r == 0:  # Middle Symmetry Boundary Condition
                        # TODO
                        radial_heatConduction = 0
                    elif r == self.radial_discretization.num_volumes:  # Reactor Wall Boundary condition
                        # TODO
                        radial_heatConduction = 0
                    else:
                        # TODO
                        radial_heatConduction = 0
                else: # 1D radial thermal conduction with U_radial = const.
                    radial_heatConduction = 4 * self.lambda_radial / self.reactorDiameter * (T[current] - self.T_wall)

                # 3.3) Combined Energy Conversation
                reactionHeat = self.reactionHeat(T[current], w_i[current, :].T, p[current])
                self.ODE_T[current] = ((
                                       - axial_convectiveHeatFlux
                                       - axial_heatConduction
                                       - radial_heatConduction
                                       #+ reactionHeat
                                       ) / (self.rho_cp_eff(w_i[current, :].T, T[current], p[current])))

            # 4) SPECIES CONSERVATION
                for comp in range(self.n_components):

                    # 4.1) Axial Species Conversation
                    if z == 0:  # Inlet Boundary Condition
                        axialMassFlow = self.axialMassFlow(T[current], w_i[current, :].T, self.w_i_in, u[current], p[current], comp) / delta_z
                    else:
                        axialMassFlow = self.axialMassFlow(T[current], w_i[current, :].T, w_i[before_z, :].T, u[current], p[current], comp) / delta_z

                    # 4.2) Radial Species Conversation
                    if self.dimension == 2:
                        if r == 0: # Middle Symmetry Boundary Condition
                            pass
                            # TODO
                        elif r == self.radial_discretization.num_volumes:  # Reactor Wall Boundary condition
                            pass
                            # TODO
                        else:
                            pass
                            # TODO

                    # 4.3) Combined Species Conversation
                    changeByReaction = self.changeByReaction(T[current], w_i[current, :].T, p[current], comp)

                    self.ODE_wi[current, comp] = ((
                                                  - axialMassFlow
                                                  - changeByReaction
                                                  ) / (self.eps * self.rho_fl(w_i[current, :].T, T[current], p[current])))


    def __reshapeCasADi_createODE(self):
        # Transform differential variables and AE/ODE in 1D vector

        w_i_reshaped = CasADi.reshape(self.w_i, self.n_spatial * self.n_components, 1)
        T_reshaped = CasADi.reshape(self.T, self.n_spatial, 1)
        x = CasADi.vertcat(w_i_reshaped, T_reshaped)

        ode_wi_reshaped = CasADi.reshape(self.ODE_wi, self.n_spatial * self.n_components, 1)
        ode_T_reshaped = CasADi.reshape(self.ODE_T, self.n_spatial, 1)
        ode = CasADi.vertcat(ode_wi_reshaped, ode_T_reshaped)

        u_reshaped = CasADi.reshape(self.u, self.n_spatial, 1)
        p_reshaped = CasADi.reshape(self.p, self.n_spatial, 1)
        z = CasADi.vertcat(u_reshaped, p_reshaped)

        ae_m_reshaped = CasADi.reshape(self.AE_m, self.n_spatial, 1)
        ae_p_reshaped = CasADi.reshape(self.AE_p, self.n_spatial , 1)
        ae = CasADi.vertcat(ae_m_reshaped, ae_p_reshaped)

        # create DAE
        self.DAE = {'x': x, 'z': z, 'ode': ode, 'alg': ae}
