from classes.ReactorSpecificQuantities.Reaction.ReactionRate import ReactionRate

class Reaction:
    def __init__(self, log):
        self.log = log
        self.reactionRate = ReactionRate(self.log)
        self.stoichiometry_coefficients = []

    def addStoichiometryCoefficient(self, name, coefficient):
        self.stoichiometry_coefficients.append([name, coefficient])

    def getReactionRate(self, T, p_CH4, p_H2O, p_CO2, p_H2):
        return self.reactionRate.rate_equation(T, p_CH4, p_H2O, p_CO2, p_H2)

    def getStoichiometryCoefficient(self, name):
        for coefficient in self.stoichiometry_coefficients:
            if name == coefficient[0]:
                return coefficient[1]
        return None
