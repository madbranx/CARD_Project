from classes.ReactorSpecificQuantities.Component.Component import Component
from classes.ReactorSpecificQuantities.Parameter.Parameter import Parameter
from classes.ReactorSpecificQuantities.Reaction.Reaction import Reaction


class ReactorSpecificQuantities:
    def __init__(self, log):
        self.log = log
        self.parameters = []
        self.components = []
        self.reaction = None

    def addParameter(self, name, value):
        parameter = Parameter(self.log, name, value)
        self.parameters.append(parameter)
        return parameter

    def addComponent(self, name):
        component = Component(self.log, name)
        self.components.append(component)
        return component

    def addReaction(self):
        reaction = Reaction(self.log)
        self.reaction = reaction
        return reaction

    def getParameterValue(self, name):
        for parameter in self.parameters:
            if parameter.getName() == name:
                return parameter.getValue()
        return None

    def getParameterSXsym(self, name):
        for parameter in self.parameters:
            if parameter.getName() == name:
                return parameter.getSXsym()
        return None

    def getComponent(self, name):
        for component in self.components:
            if component.getName() == name:
                return component
        return None

    def getComponents(self):
        return self.components

    def getReactionRate(self):
        return self.reaction.getReactionRate()

    def getStoichCoeffs(self):
        return self.reaction.getStoichiometricCoefficients()



