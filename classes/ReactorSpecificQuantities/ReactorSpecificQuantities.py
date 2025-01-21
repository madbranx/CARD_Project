from classes.ReactorSpecificQuantities.Component.Component import Component
from classes.ReactorSpecificQuantities.Parameter.Parameter import Parameter
from classes.ReactorSpecificQuantities.Reaction.Reaction import Reaction


class ReactorSpecificQuantities:
    def __init__(self, log):
        self.log = log
        self.parameters = []
        self.components = []
        self.catalyst = None
        self.reaction = None

    def addParameter(self, name, value):
        self.log.addEntry("adding Parameter " + name + " = " + str(value), 2)
        parameter = Parameter(self.log, name, value)
        self.parameters.append(parameter)
        return parameter

    def addComponent(self, name):
        self.log.addEntry("adding Component " + name, 2)
        component = Component(self.log, name)
        self.components.append(component)
        return component

    def addCatalyst(self, name):
        self.log.addEntry("adding Catalyst " + name, 2)
        catalyst = Component(self.log, name)
        self.catalyst = catalyst
        return catalyst

    def addReaction(self):
        self.log.addEntry("adding Reaction", 2)
        reaction = Reaction(self.log)
        self.reaction = reaction
        return reaction

    def calculate_void_fraction(self):
        d_cat = self.getParameterValue("cat_diameter")
        d_reactor = self.getParameterValue("reactorDiameter")

        ratio = d_cat / d_reactor

        if ratio <= 0.5:
            epsilon = 0.4 + 0.05 * ratio + 0.412 * (ratio ** 2)
        elif 0.5 < ratio <= 0.536:
            epsilon = 0.528 + 2.464 * (ratio - 0.5)
        else:  # ratio >= 0.536
            epsilon = 1 - 0.667 * (ratio ** 3) * (2 * ratio - 1) ** -0.5

        return epsilon

    def getParameterValue(self, name):
        for parameter in self.parameters:
            if parameter.getName() == name:
                return parameter.getValue()
        return None

    def getComponent(self, name):
        for component in self.components:
            if component.getName() == name:
                return component
        return None

    def getComponents(self):
        return self.components

    def getNComponents(self):
        return len(self.components)

    def getCatalyst(self):
        return self.catalyst

    def getReactionRate(self):
        return self.reaction.getReactionRate()

    def addStoichCoeff(self, name, coefficient):
        self.reaction.addStoichiometryCoefficient(name, coefficient)

    def getStoichCoeff(self, name):
        return self.reaction.getStoichiometryCoefficient(name)

    def getStoichCoeffs(self):
        return self.reaction.getStoichiometryCoefficients()

    def getMolarWeights(self):
        Mw_i = []
        for component in self.getComponents():
            Mw_i.append(component.get_molecular_weight())
        return Mw_i


