from .EnergyConservation.EnergyConservation import EnergyConservation
from .MassConservation.MassConservation import MassConservation
from .PressureDrop.PressureDrop import PressureDrop
from .SpeciesConservation.SpeciesConservation import SpeciesConservation
from ..Discretization.Discretization import Discretization
from ..ReactorSpecificQuantities.ReactorSpecificQuantities import ReactorSpecificQuantities
from ..ReactorSpecificQuantities.Component.Component import Component, MaterialProperty

class FixedBedReactor:
    ONE_D = 1
    TWO_D = 2
    def __init__(self, log, dimension, n_z, n_r=None):
        self.log = log
        # Discretisation
        self.dimension = dimension
        self.disc_z = None
        self.n_z = n_z
        if dimension == self.TWO_D:
            self.disc_r = None
            self.n_r = n_r
        # Reactor Specific Quantities
        self.RSQ = None
        # Conservations
        self.SpeciesConservation = None
        self.EnergyConservation = None
        self.MassConservation = None
        self.PressureDrop = None
        # SetUp
        self.__setUp()

    def __setUp(self):
        self.__initiateRSQ()
        self.__create_spacialDiscretization()
        self.__create_conservationEquations()
        self.__create_DAEstruct()

    def __initiateRSQ(self):
        self.RSQ = ReactorSpecificQuantities(self.log)
        # add Parameters to RSQ
        self.RSQ.addParameter("reactorLength", 10)
        self.RSQ.addParameter("reactorDiameter", 1)

        # add Components to RSQ
        H2 = self.RSQ.addComponent("H2")
        H2.add_property(Component.DENSITY, 100)
        H2.add_property(Component.THERMAL_CONDUCTIVITY, [1, 0.1], MaterialProperty.LINEAR)

        # add Reaction to RSQ
        self.RSQ.addReaction()
        
    def __create_spacialDiscretization(self):
        reactorLength = self.RSQ.getParameterValue("reactorLength")
        reactorRadius = self.RSQ.getParameterValue("reactorDiameter") / 2
        self.disc_z = Discretization(self.log, self.n_z, Discretization.EQUIDISTANT, start=0.0, end=reactorLength)
        if self.dimension == self.TWO_D:
            self.disc_r = Discretization(self.log, self.n_r, Discretization.EQUIDISTANT, start=0.0, end=reactorRadius)

    def __create_conservationEquations(self):
        self.speciesConservation = SpeciesConservation(self.log)
        self.EnergyConservation = EnergyConservation(self.log)
        self.MassConservation = MassConservation(self.log)
        self.PressureDrop = PressureDrop(self.log)

    def __create_DAEstruct(self):
        pass  # TODO
