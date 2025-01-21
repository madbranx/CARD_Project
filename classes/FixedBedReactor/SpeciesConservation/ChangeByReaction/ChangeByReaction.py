from classes.FixedBedReactor.SpeciesConservation.ChangeByReaction.EffectivenessFactor.EffectivenessFactor import \
    EffectivenessFactor


class ChangeByReaction:

    def __init__(self,log, RSQ, GCF):
        self.log = log
        self.RSQ = RSQ
        self.GCF = GCF
        self.effFactor = EffectivenessFactor(self.log, self.RSQ, self.GCF)

    # Following Methods use location specific Arguments -> can be used for 1D and 2D!
    def calc(self, T, w_i, comp):
        eps = self.RSQ.getParameterValue("bed_void_fraction")
        Mw_i = self.RSQ.getMolarWeights()
        stoich_coeffs = self.RSQ.getStoichCoeffs()

        #TODO use get Reaction -> Reaction.calcReactionrate
        # calcReactionrate should use w_i as input -> define function in GCF to calc p from w_i
        # give reaction rate to eff factor as argument
        # .
        # Dont use RSQ to get RRate/Stoich ....
        # Use Reaction directly -> iterate over reactions: for reaction in self.RSQ.getReactions(): ....

        reaction_rate = self.RSQ.getReactionRate()
        eff_factor = self.effFactor.calc(T, w_i)

        return (1-eps)*Mw_i[comp]*stoich_coeffs[comp]*eff_factor*reaction_rate

