from classes.ReactorSpecificQuantities.Reaction.ReactionRateKoschany import ReactionRateKoschany

class Reaction:
    def __init__(self, log):
        self.log = log
        self.reactionRate = None
        self.stoichiometry_coefficients = None

    def setReactionRate(self, reactionRate):
        self.reactionRate = reactionRate

    def setStoichiometryCoefficients(self, coefficients):
        self.stoichiometry_coefficients = coefficients

    def getStoichiometryCoefficients(self):
        return self.stoichiometry_coefficients

    def getReactionRate(self, w_i, T):
        return self.reactionRate.rate_equation(w_i, T)
