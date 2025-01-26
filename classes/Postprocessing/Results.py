class Results:
    def __init__(self, axialDiscretization, radialDiscretization, timeDiscretization):
        self.axialDiscretization = axialDiscretization
        self.radialDiscretization = radialDiscretization
        self.timeDiscretization = timeDiscretization

        self.w_i = None
        self.T = None
        self.p = None
        self.u = None

    def add_values(self, w_i, T, p, u):
        self.w_i = w_i
        self.T = T
        self.p = p
        self.u = u

    def get_values(self, t_step):
        return self.w_i[:, t_step, :], self.T[:, t_step], self.p[:, t_step], self.u[:, t_step]

    def get_z_pos(self):
        z_pos = self.axialDiscretization.get_centroids()
        return z_pos




