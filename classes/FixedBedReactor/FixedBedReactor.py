from ..ReactorSpecificQuantities.ReactorSpecificQuantities import ReactorSpecificQuantities
from ..ReactorSpecificQuantities.Component.Component import Component, MaterialProperty

class FixedBedReactor:
    def __init__(self, log):
        self.log = log
        self.RSQ = ReactorSpecificQuantities(self.log)
        self.__initiateRSQ()

    def __initiateRSQ(self):
        # add Parameters to RSQ
        self.RSQ.addParameter("reactorLength", 1)

        # add Components to RSQ
        H2 = self.RSQ.addComponent("H2")
        H2.add_property(Component.DENSITY, 100)
        H2.add_property(Component.THERMAL_CONDUCTIVITY, [1, 0.1], MaterialProperty.LINEAR)

        # add Reaction to RSQ
        self.RSQ.addReaction()


