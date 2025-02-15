import numpy as np

class Results:
    def __init__(self, axialDiscretization, radialDiscretization, timeDiscretization):
        self.axialDiscretization = axialDiscretization
        self.radialDiscretization = radialDiscretization
        self.timeDiscretization = timeDiscretization

        self.w_i = None
        self.T = None
        self.p = None
        self.u = None

        self.reactor = None

    def set_reactor(self, reactor):
        self.reactor = reactor

    def add_values(self, w_i, T, u, p):
        self.w_i = w_i
        self.T = T
        self.p = p
        self.u = u

    def get_values(self, t_step):
        return self.w_i[:, t_step, :], self.T[:, t_step], self.p[:, t_step], self.u[:, t_step]

    def get_z_pos(self):
        z_pos = self.axialDiscretization.get_centroids()
        return z_pos

    def get_r_pos(self):
        r_pos = self.radialDiscretization.get_centroids()
        return r_pos

    def get_r_faces(self):
        r_faces = self.radialDiscretization.get_faces()
        return r_faces

    def get_2D_values(self, t_step):
        # TODO number of components hardcoded
        n_axial_volumes = self.axialDiscretization.num_volumes
        n_radial_volumes = self.radialDiscretization.num_volumes

        w_i_2D = np.zeros((4, n_axial_volumes, n_radial_volumes))
        T_2D = np.zeros((n_axial_volumes, n_radial_volumes))
        p_2D = np.zeros((n_axial_volumes, n_radial_volumes))
        u_2D = np.zeros((n_axial_volumes, n_radial_volumes))

        # Convert results vector into axial and radial accessible matrix
        for r in range(n_radial_volumes):
            for z in range(n_axial_volumes):
                current = z + r * n_axial_volumes

                w_i_2D[:,z,r] = self.w_i[current,t_step,:]
                T_2D[z,r] = self.T[current,t_step]
                p_2D[z,r] = self.p[current,t_step]
                u_2D[z,r] = self.u[current,t_step]

        return w_i_2D, T_2D, p_2D, u_2D

    def getConversion_2D(self, t_step, comp):
        n_axial_volumes = self.axialDiscretization.num_volumes
        n_radial_volumes = self.radialDiscretization.num_volumes

        i = comp
        w_i_in = self.reactor.w_i_in
        w_M = self.reactor.getMolarWeights()
        tot_moles = 0
        for i, mw in enumerate(w_M):
            tot_moles += w_i_in[i]/mw

        n_i_in = w_i_in[i]/w_M[i]/tot_moles

        X_2D = np.zeros((n_axial_volumes, n_radial_volumes))

        for r in range(n_radial_volumes):
            for z in range(n_axial_volumes):
                current = z + r * n_axial_volumes
                w_i =  self.w_i[current, t_step, :]
                n_i = self.reactor.moleFractions(w_i)
                X_2D[z, r] = (n_i_in-n_i[i])/n_i_in

        return X_2D

    def average_trapezoidal(self, f_values):
        r_faces = self.radialDiscretization.get_faces()
        R = r_faces[-1]  # Maximum radius
        N = len(f_values)  # Number of cells

        integral = 0.0
        for i in range(N):
            r_mid = (r_faces[i] + r_faces[i + 1]) / 2  # Midpoint radius
            dr = r_faces[i + 1] - r_faces[i]  # Cell width
            integral += f_values[i] * r_mid * dr  # Trapezoidal rule

        avg_f = (2 * np.pi * integral) / (np.pi * R ** 2)  # Normalize by cross-section area
        return avg_f


    def get_rawValues(self, timestep):
        return self.w_i[:, timestep, :], self.T[:, timestep], self.p[:, timestep], self.u[:, timestep]


    def getConversion_1D(self, t_step, comp):
        n_axial_volumes = self.axialDiscretization.num_volumes

        i = comp
        w_i_in = self.reactor.w_i_in
        w_M = self.reactor.getMolarWeights()
        tot_moles = 0
        for i, mw in enumerate(w_M):
            tot_moles += w_i_in[i]/mw

        n_i_in = w_i_in[i]/w_M[i]/tot_moles

        X = np.zeros(n_axial_volumes)

        for z in range(n_axial_volumes):
            w_i =  self.w_i[z, t_step, :]
            n_i = self.reactor.moleFractions(w_i)
            X[z] = (n_i_in-n_i[i])/n_i_in

        return X


