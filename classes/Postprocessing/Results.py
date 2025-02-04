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




