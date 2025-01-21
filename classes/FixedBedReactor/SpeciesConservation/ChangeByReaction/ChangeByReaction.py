from classes.FixedBedReactor.SpeciesConservation.ChangeByReaction.EffectivenessFactor.EffectivenessFactor import \
    EffectivenessFactor

class ChangeByReaction:

    def __init__(self,log, RSQ, GCF):
        self.log = log
        self.RSQ = RSQ
        self.GCF = GCF
        self.effFactor = EffectivenessFactor(self.log, self.RSQ, self.GCF)

    # Following Methods use location specific Arguments -> can be used for 1D and 2D!
    def calc(self, T, w_i, p, comp):
        eps = self.RSQ.getParameterValue("bed_void_fraction")
        Mw_i = self.RSQ.getMolarWeights()

        # effectiveness Factor is always calculated with the reaction 1 -> no function of the other reactions
        eff_factor = self.effFactor.calc(w_i, T, p)

        reactions = self.RSQ.getReactions()
        sum_stoichiometric_rate = 0

        # calculate the stoichiometric factor & reaction rate for each reaction and sum them up
        for reaction in reactions:
            stoich_factor = reaction.getStoichiometryCoefficients()[comp]
            reaction_rate = reaction.getReactionRate(w_i, T)
            sum_stoichiometric_rate = (sum_stoichiometric_rate + stoich_factor * reaction_rate)

        return (1-eps)*Mw_i[comp]*sum_stoichiometric_rate*eff_factor

