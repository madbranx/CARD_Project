from ReactionRate import *

class Reaction:
    def __init__(self, log):
        self.log = log
        self.reactionRate = ReactionRate(self.log)



