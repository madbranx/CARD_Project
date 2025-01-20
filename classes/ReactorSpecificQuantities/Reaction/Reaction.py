from classes.ReactorSpecificQuantities.Reaction.ReactionRate import ReactionRate

class Reaction:
    def __init__(self, log):
        self.log = log
        self.reactionRate = ReactionRate(self.log)
        self.stoichiometry_coefficients = []

    def setStoichiometryCoefficients(self, stoichiometry_coefficients):
        self.stoichiometry_coefficients = stoichiometry_coefficients

    def getReactionRate(self, T, p_CH4, p_H2O, p_CO2, p_H2):
        return self.reactionRate.rate_equation(T, p_CH4, p_H2O, p_CO2, p_H2)

    def getStoichiometryCoefficients(self):
        return self.stoichiometry_coefficients
