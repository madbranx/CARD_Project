class Reaction:
    def __init__(self, log):
        self.log = log
        self.reactionRate = None
        self.stoichiometry_coefficients = None
        self.reaction_Enthalpy = 0

    def setReactionRate(self, reactionRate):
        self.reactionRate = reactionRate

    def setStoichiometryCoefficients(self, coefficients):
        self.stoichiometry_coefficients = coefficients

    def setReactionEnthalpy(self, enthalpy):
        self.reaction_Enthalpy = enthalpy

    def getStoichiometryCoefficients(self):
        return self.stoichiometry_coefficients

    def getReactionRate(self, w_i, T, p):
        return self.reactionRate.rate_equation(w_i, T, p)

    def getReactionEnthalpy(self):
        return self.reaction_Enthalpy
