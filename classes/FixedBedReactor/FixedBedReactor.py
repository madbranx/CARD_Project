from .EnergyConservation.EnergyConservation import EnergyConservation
from .MassConservation.MassConservation import MassConservation
from .PressureDrop.PressureDrop import PressureDrop
from .SpeciesConservation.MolecularDiffCoff import MolecularDiffCoff
from .SpeciesConservation.SpeciesConservation import SpeciesConservation
from ..Discretization.Discretization import Discretization
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
        if dimension == self.TWO_D:
            self.disc_r = None
            self.n_radial = n_radial
        # Reactor Specific Quantities
        self.RSQ = None
        # Conservations
        self.SpeciesConservation = None
        self.EnergyConservation = None
        self.MassConservation = None
        self.PressureDrop = None

    def setUp(self):
        self.__initiateRSQ()
        self.__create_spacialDiscretization()
        self.__create_conservationEquations()
        self.__create_DAEstruct()
        self.log.updateLog()

    def __initiateRSQ(self):
        self.log.addEntry("Initiating RSQ", 1)
        self.RSQ = ReactorSpecificQuantities(self.log)

        # add Parameters to RSQ
        self.RSQ.addParameter("reactorLength", 10)
        self.RSQ.addParameter("reactorDiameter", 1)

        # add Components to RSQ
        H2 = self.RSQ.addComponent("H2")
        H2.add_property(Component.DENSITY, 100)
        H2.add_property(Component.THERMAL_CONDUCTIVITY, [1, 0.1], MaterialProperty.LINEAR)
        H2.add_property(Component.HEAT_CAPACITY, [1, 0.1, 0.01], MaterialProperty.POLYNOMIAL)

        #add Catalyst to RSQ
        cat = self.RSQ.addCatalyst("cat1")
        cat.add_property(Component.DENSITY, 2000)

        # add Reaction to RSQ
        self.RSQ.addReaction()
        
    def __create_spacialDiscretization(self):
        self.log.addEntry("Creating spatial discretization for " + str(self.dimension) + " dimension(s)", 1)

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
        self.log.addEntry("Creating conservation equations", 1)

        MDC = MolecularDiffCoff() #TODO

        self.speciesConservation = SpeciesConservation(self.log, self.dimension, self.RSQ)
        self.EnergyConservation = EnergyConservation(self.log, self.dimension, self.RSQ)
        self.MassConservation = MassConservation(self.log, self.dimension, self.RSQ)
        self.PressureDrop = PressureDrop(self.log, self.dimension, self.RSQ)

    def __create_DAEstruct(self):
        # TODO differentiate between 1/2D

        n_comp = self.RSQ.getNComponents() # get from RSQ
        print(n_comp)

        if self.dimension == FixedBedReactor.ONE_D:
            w_i = CasADi.SX.sym('w_i', (self.n_axial, n_comp))  # A,B mol/m^3
            T = CasADi.SX.sym('T', self.n_axial)  # [K]
            # Define algebraic variables
            u = CasADi.SX.sym('u', self.n_axial)  # [m/s]
        elif self.dimension == FixedBedReactor.TWO_D:
            w_i = CasADi.SX.sym('w_i', (self.n_axial, self.n_radial, n_comp))  # A,B mol/m^3
            T = CasADi.SX.sym('T', self.n_axial, self.n_radial)  # [K]
            # Define algebraic variables
            u = CasADi.SX.sym('u', self.n_axial,self.n_radial)  # [m/s]
        else:
            self.log.addError("dimension not supported")
            return None

        # create AEs/ODEs
        ae_mass = CasADi.SX.sym("ae_mass", self.n_axial)
        self.MassConservation.createCasADi(ae_mass, T, w_i, u)
