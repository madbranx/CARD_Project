from .EnergyConservation.EnergyConservation import EnergyConservation
from .MassConservation.MassConservation import MassConservation
from .PressureDrop.PressureDrop import PressureDrop
from .SpeciesConservation.SpeciesConservation import SpeciesConservation
from ..Discretization.Discretization import Discretization
from ..GeneralConversionFunctions.GeneralConversionFunctions import GeneralConversionFunctions
from ..ReactorSpecificQuantities.Reaction.Reaction import Reaction
from ..ReactorSpecificQuantities.Reaction.ReactionRateKoschany import ReactionRateKoschany
from ..ReactorSpecificQuantities.ReactorSpecificQuantities import ReactorSpecificQuantities
from ..ReactorSpecificQuantities.Component.Component import Component, MaterialProperty

import casadi as CasADi

class FixedBedReactor:
    ONE_D = 1
    TWO_D = 2
    def __init__(self, log, dimension, n_axial, n_radial=None):
        self.log = log
        log.addEntry("creating FixedBedReactor instance", 0)

        # Discretization
        self.dimension = dimension
        self.disc_z = None
        self.n_axial = n_axial
        self.disc_r = None
        self.n_radial = n_radial

        # Reactor Specific Quantities
        self.RSQ = None
        self.n_components = None

        self.GCF = None

        # Conservations
        self.SpeciesConservation = None
        self.EnergyConservation = None
        self.MassConservation = None
        self.PressureDrop = None
        self.dae = None

    def setUp(self):
        self.__initiate_RSQ_GCF()
        self.__create_spacialDiscretization()
        self.__create_conservationEquations()
        self.__create_DAEstruct()
        self.log.updateLog()

    def __initiate_RSQ_GCF(self):
        self.log.addEntry("Initiating RSQ", 1)
        self.RSQ = ReactorSpecificQuantities(self.log)
        self.GCF = GeneralConversionFunctions(self.log, self.RSQ)

        # TODO add remaining Parameters/Properties

        # add Constants to RSQ
        self.RSQ.addParameter("R", 8.314)

        # add Parameters to RSQ
        self.RSQ.addParameter("reactorLength", 10)
        self.RSQ.addParameter("reactorDiameter", 1)

        self.RSQ.addParameter("cat_diameter", 1)
        self.RSQ.addParameter("cat_tortuosity", 1)
        self.RSQ.addParameter("cat_porosity", 1)

        self.RSQ.addParameter("diameter_pore", 1)

        # calculate void fraction
        eps = self.RSQ.calculate_void_fraction()
        self.RSQ.addParameter("bed_void_fraction", eps)

        self.RSQ.addParameter("T_in", 200)
        self.RSQ.addParameter("u_in", 1)
        self.RSQ.addParameter("w_i_in", [0, 1, 0, 0])
        self.RSQ.addParameter("p_in", 1e5)

        # add Components and their Properties to RSQ
        CH4 = self.RSQ.addComponent("CH4")
        CH4.add_property(Component.DENSITY, 1)
        CH4.add_property(Component.HEAT_CAPACITY, 1)
        CH4.add_property(Component.THERMAL_CONDUCTIVITY, 1)
        CH4.add_property(Component.COLLISION_AREA, 1)
        CH4.add_property(Component.DIFFUSION_VOLUME, 1)
        CH4.add_property(Component.DYNAMIC_VISCOSITY, 1)
        CH4.add_property(Component.MOLECULAR_WEIGHT, 1)

        H20 = self.RSQ.addComponent("H2O")
        H20.add_property(Component.DENSITY, 2)
        H20.add_property(Component.HEAT_CAPACITY, 2)
        H20.add_property(Component.THERMAL_CONDUCTIVITY, 2)
        H20.add_property(Component.COLLISION_AREA, 2)
        H20.add_property(Component.DIFFUSION_VOLUME, 2)
        H20.add_property(Component.DYNAMIC_VISCOSITY, 2)
        H20.add_property(Component.MOLECULAR_WEIGHT, 2)

        CO2 = self.RSQ.addComponent("CO2")
        CO2.add_property(Component.DENSITY, 3)
        CO2.add_property(Component.HEAT_CAPACITY, 3)
        CO2.add_property(Component.THERMAL_CONDUCTIVITY, 3)
        CO2.add_property(Component.COLLISION_AREA, 3)
        CO2.add_property(Component.DIFFUSION_VOLUME, 3)
        CO2.add_property(Component.DYNAMIC_VISCOSITY, 3)
        CO2.add_property(Component.MOLECULAR_WEIGHT, 3)

        H2 = self.RSQ.addComponent("H2")
        H2.add_property(Component.DENSITY, 4)
        H2.add_property(Component.HEAT_CAPACITY, 4)
        H2.add_property(Component.THERMAL_CONDUCTIVITY, 4)
        H2.add_property(Component.COLLISION_AREA, 4)
        H2.add_property(Component.DIFFUSION_VOLUME, 4)
        H2.add_property(Component.DYNAMIC_VISCOSITY, 4)
        H2.add_property(Component.MOLECULAR_WEIGHT, 4)

        #add Catalyst to RSQ
        cat = self.RSQ.addCatalyst("cat1")
        cat.add_property(Component.DENSITY, 2000)

        # add Reactions to RSQ

        # CH4 + 2 H20 -> CO2 + 4 H2
        reaction1 = Reaction(self.log)
        reaction1.setStoichiometryCoefficients([-1, -2, 1, 4])
        # set reaction Rate
        reactionRate1 = ReactionRateKoschany(self.log, self.GCF, self.RSQ)
        reaction1.setReactionRate(reactionRate1)

        self.RSQ.addReaction(reaction1)

        # add component number to class attribute
        self.n_components = self.RSQ.getNComponents()  # get from RSQ

    def __create_spacialDiscretization(self):
        self.log.addEntry("creating spatial discretization for " + str(self.dimension) + " dimension(s)", 1)

        # Axial discretization
        reactorLength = self.RSQ.getParameterValue("reactorLength")
        self.log.addEntry("axial discretization", 2)
        #self.disc_z = Discretization(self.log, self.n_z, Discretization.EQUIDISTANT, start=0.0, end=reactorLength)
        ranges = [(0, reactorLength/2, 4), (reactorLength/2, reactorLength, 1)]
        self.disc_z = Discretization(self.log, self.n_axial, Discretization.ARRAY, start=0.0, end=reactorLength, ranges=ranges)

        # Radial discretization
        reactorRadius = self.RSQ.getParameterValue("reactorDiameter") / 2
        if self.dimension == self.TWO_D:
            self.log.addEntry("radial  discretization", 2)
            self.disc_r = Discretization(self.log, self.n_radial, Discretization.EQUIDISTANT, start=0.0, end=reactorRadius)

    def __create_conservationEquations(self):
        self.log.addEntry("creating conservation equations", 1)

        self.SpeciesConservation = SpeciesConservation(self.log, self.dimension, self.RSQ, self.GCF, self.disc_z, self.disc_r)
        self.EnergyConservation = EnergyConservation(self.log, self.dimension, self.RSQ, self.GCF)
        self.MassConservation = MassConservation(self.log, self.dimension, self.RSQ, self.GCF)
        self.PressureDrop = PressureDrop(self.log, self.dimension, self.RSQ, self.GCF, self.disc_z, self.disc_r)

    def __create_DAEstruct(self):
        self.log.addEntry("creating CasADi DAE structure", 1)
        [w_i, T, u, p] = self.__initialize_CasADi_symbols()

        # create AEs/ODEs
        self.log.addEntry("initializing CasADi AE and ODE structures", 2)
        if self.dimension == FixedBedReactor.ONE_D:
            ae_m= CasADi.SX.sym("ae_mass", self.n_axial)
            ae_p = CasADi.SX.sym("ae_pressure_drop", self.n_axial)
            ode_wi = CasADi.SX.sym("ode_wi", (self.n_axial, self.n_components))
            ode_T = CasADi.SX.sym("ode_wi", self.n_axial)

        elif self.dimension == FixedBedReactor.TWO_D:
            ae_m = CasADi.SX.sym("ae_mass", self.n_axial* self.n_radial)
            ae_p = CasADi.SX.sym("ae_pressure_drop", self.n_axial* self.n_radial)
            ode_wi = CasADi.SX.sym("ode_wi", (self.n_axial* self.n_radial, self.n_components))
            ode_T = CasADi.SX.sym("ode_wi", self.n_axial* self.n_radial)

        else:
           return None

        # Define AEs/ODEs
        self.log.addEntry("defining CasADi AE and ODE equations", 2)
        self.MassConservation.createCasADi(ae_m, T, w_i, u, p)
        self.PressureDrop.createCasADi(ae_p, T, w_i, u, p)
        self.SpeciesConservation.createCasADi(ode_wi, T, w_i, u, p)
        self.EnergyConservation.createCasADi(ode_T, T, w_i, u, p)

        # Reshape AEs/ODEs
        self.log.addEntry("reshaping CasADi AE and ODE structures", 2)
        [ae, ode, x, z] = self.__reshape_CasADi_ae_ode(w_i, T, u, p, ae_m, ae_p, ode_wi, ode_T)

        # create DAE struct
        self.log.addEntry("creating CasADi DAE object", 2)
        self.dae = {'x': x, 'z': z, 'ode': ode, 'alg': ae}


    def __initialize_CasADi_symbols(self):
        self.log.addEntry("initializing CasADi symbols", 2)
        if self.dimension == FixedBedReactor.ONE_D:
            # Define differential variables
            w_i = CasADi.SX.sym('w_i', (self.n_axial, self.n_components))
            self.log.addEntry("w_i = " + str(self.n_axial) + " x " + str(self.n_components), 3)
            T = CasADi.SX.sym('T', self.n_axial)
            self.log.addEntry("T = " + str(self.n_axial), 3)

            # Define algebraic variables
            u = CasADi.SX.sym('u', self.n_axial)
            self.log.addEntry("u = " + str(self.n_axial), 3)
            p = CasADi.SX.sym('p', self.n_axial)
            self.log.addEntry("p = " + str(self.n_axial), 3)

        elif self.dimension == FixedBedReactor.TWO_D:
            # Casadi does not support 3D elements -> axial and radial elements are combined (self.n_axial* self.n_radial)!

            # Define differential variables
            w_i = CasADi.SX.sym('w_i', (self.n_axial* self.n_radial, self.n_components))
            self.log.addEntry("w_i = (" + str(self.n_axial) + " * " + str(self.n_radial) + ") x " + str(self.n_components), 3)
            T = CasADi.SX.sym('T', (self.n_axial* self.n_radial))
            self.log.addEntry("T = (" + str(self.n_axial) + " * " + str(self.n_radial) + ")", 3)

            # Define algebraic variables
            u = CasADi.SX.sym('u', (self.n_axial* self.n_radial))
            self.log.addEntry("u = (" + str(self.n_axial) + " * " + str(self.n_radial) + ")", 3)
            p = CasADi.SX.sym('p', (self.n_axial* self.n_radial))
            self.log.addEntry("p = (" + str(self.n_axial) + " * " + str(self.n_radial) + ")", 3)

        else:
            self.log.addError("dimension not supported")
            return None

        return w_i, T, u, p


    def __reshape_CasADi_ae_ode(self,w_i, T, u, p, ae_m, ae_p, ode_wi, ode_T):
        # Transform differential variable and ode in 1D vector

        if self.dimension == FixedBedReactor.ONE_D:
            w_i_reshaped = CasADi.reshape(w_i, self.n_axial * self.n_components, 1)
            x = CasADi.vertcat(w_i_reshaped, T)

            ode_wi_reshaped = CasADi.reshape(ode_wi, self.n_axial * self.n_components, 1)
            ode = CasADi.vertcat(ode_wi_reshaped, ode_T)

            ae = CasADi.vertcat(ae_m, ae_p)
            z = CasADi.vertcat(u, p)

        elif self.dimension == FixedBedReactor.TWO_D:
            w_i_reshaped = CasADi.reshape(w_i, self.n_axial * self.n_components * self.n_radial, 1)
            T_reshaped = CasADi.reshape(T, self.n_axial * self.n_radial, 1)
            x = CasADi.vertcat(w_i_reshaped, T_reshaped)

            ode_wi_reshaped = CasADi.reshape(ode_wi, self.n_axial * self.n_components * self.n_radial, 1)
            ode_T_reshaped = CasADi.reshape(ode_T, self.n_axial * self.n_radial, 1)
            ode = CasADi.vertcat(ode_wi_reshaped, ode_T_reshaped)

            u_reshaped = CasADi.reshape(u, self.n_axial * self.n_radial, 1)
            p_reshaped = CasADi.reshape(p, self.n_axial * self.n_radial, 1)
            z = CasADi.vertcat(u_reshaped, p_reshaped)

            ae_m_reshaped = CasADi.reshape(ae_m, self.n_axial * self.n_radial, 1)
            ae_p_reshaped = CasADi.reshape(ae_p, self.n_axial * self.n_radial, 1)
            ae = CasADi.vertcat(ae_m_reshaped, ae_p_reshaped)

        else:
            self.log.addError("dimension not supported")
            return None

        return ae, ode, x, z
