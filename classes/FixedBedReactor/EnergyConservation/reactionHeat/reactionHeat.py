from classes.FixedBedReactor.SpeciesConservation.ChangeByReaction.EffectivenessFactor.EffectivenessFactor import \
    EffectivenessFactor
from classes.ReactorSpecificQuantities.Component.Component import Component


class ReactionHeat:

    def __init__(self,log, GCF, RSQ):
        self.log = log
        self.GCF = GCF
        self.RSQ = RSQ
        self.effFactor = EffectivenessFactor(self.log, self.RSQ, self.GCF)

    # Following Methods use location specific Arguments -> can be used for 1D and 2D!
    def calc(self, T, w_i, p):
        eps = self.RSQ.getParameterValue("bed_void_fraction")

        # effectiveness Factor is always calculated with the reaction 1 -> no function of the other reactions
        eff_factor = self.effFactor.calc(w_i, T, p)

        # calculate the stoichiometric factor & reaction rate for each reaction and sum them up
        sum_rate_deltaH = 0
        reactions = self.RSQ.getReactions()
        for reaction in reactions:
            reaction_rate = reaction.getReactionRate(w_i, T, p)
            delta_H_r = reaction.getReactionEnthalpy()
            sum_rate_deltaH += delta_H_r * reaction_rate

        return (1-eps)*eff_factor*sum_rate_deltaH
