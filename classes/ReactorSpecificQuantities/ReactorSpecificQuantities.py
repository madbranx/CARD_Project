from Parameter.Parameter import Parameter
from Component.Component import Component
from Reaction.Reaction import Reaction


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

    def getParameters(self):
        return self.parameters

    def getComponents(self):
        return self.components

    def getReaction(self):
        return self.reaction

    def getParameter(self, name):
        for parameter in self.parameters:
            if parameter.getName() == name:
                return parameter
        return None

    def getComponent(self, name):
        for component in self.components:
            if component.getName() == name:
                return component
        return None



