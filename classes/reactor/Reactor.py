class Reactor:
    def __init__(self, d_reactor, l_reactor, d_particle):
        self.diameter = d_reactor
        self.length = l_reactor
        self.voidFraction = self.__calculate_voidfraction(d_reactor, l_reactor, d_particle)

    def __calculate_voidfraction(self, d_reactor, l_reactor, d_particle):
        eps = 1
        return eps

    def get_length(self):
        return self.length

    def get_diameter(self):
        return self.diameter

    def get_voidfraction(self):
        return self.voidFraction