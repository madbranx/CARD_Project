from classes.Parameters.Discretization import Discretization
from classes.Postprocessing.Results import Results

import casadi as CasADi
import numpy as np


class Integrator:
    def __init__(self, reactor):
        self.reactor = reactor

        self.time_discretization = None
        self.axial_discretization = reactor.axial_discretization
        self.radial_discretization = reactor.radial_discretization

        self.integrator = None
        self.x_0 = None
        self.z_0 = None
        self.results = None

        self.w_i_res = None
        self.T_res = None
        self.u_res = None
        self.p_res = None

    def setup(self, abstol, reltol, t_start, t_stop, t_steps):

        self.time_discretization = Discretization(t_steps, start=t_start, end=t_stop)
        options = {'abstol': abstol, 'reltol': reltol}
        dae = self.reactor.DAE
        timepoints = self.time_discretization.get_faces()

        self.integrator = CasADi.integrator('I', 'idas', dae, timepoints[0], timepoints, options)
        self.__setInitialValues()

    def __setInitialValues(self):
        w_i_in = self.reactor.w_i_in
        T_in = self.reactor.T_in
        u_in = self.reactor.u_in
        p_in = self.reactor.p_in

        n_axial = self.axial_discretization.num_volumes
        n_radial = self.radial_discretization.num_volumes
        n_components = self.reactor.n_components

        x_0 = np.zeros(n_axial * n_radial * (n_components + 1))

        # setting initial values for w_i
        for comp in range(n_components):
            if comp == 0:
                x_0[0: n_axial*n_radial] = w_i_in[0]
            else:
                x_0[n_axial*n_radial*comp : n_axial*n_radial*(comp+1)] = w_i_in[comp]

        # setting initial values for T
        x_0[n_components*n_axial*n_radial:] = T_in
        self.x_0 = x_0

        z_0= np.zeros(n_axial * n_radial * 2)
        z_0[0:n_axial * n_radial] = u_in
        z_0[n_axial * n_radial:] = p_in
        self.z_0 = z_0

    def integrate(self):
        print("starting integration ...")
        self.results =  self.integrator(x0=self.x_0, z0=self.z_0)
        self.__extractResults()
        self.__printMassDeviation()

        results = Results(self.axial_discretization, self.radial_discretization, self.time_discretization)
        results.add_values(self.w_i_res, self.T_res, self.u_res, self.p_res)

        return results

    def __extractResults(self):
        n_spatial = self.axial_discretization.num_volumes * self.radial_discretization.num_volumes
        n_components = self.reactor.n_components
        t_steps = self.time_discretization.num_volumes+1

        w_i_res = np.empty(shape=(n_spatial, t_steps, n_components))
        T_res = np.empty(shape=(n_spatial, t_steps))
        res_x = self.results['xf'].full()

        for t in range(t_steps):
            T_res[:, t] = res_x[n_components * n_spatial:, t]
            for comp in range(n_components):
                w_i_res[:, t, comp] = res_x[comp * n_spatial: n_spatial * (comp + 1), t]

        self.w_i_res = w_i_res
        self.T_res = T_res

        ae_res = self.results['zf'].full()

        self.u_res = ae_res[:n_spatial, :]
        self.p_res = ae_res[n_spatial:, :]

    def __printMassDeviation(self):
        n_spatial = self.axial_discretization.num_volumes * self.radial_discretization.num_volumes
        t_steps = self.time_discretization.num_volumes+1

        MassFluxDev = np.empty(shape=(n_spatial, t_steps))
        mdot_0 = (self.reactor.u_in * self.reactor.rho_fl(self.reactor.w_i_in, self.reactor.T_in, self.reactor.p_in)).__float__()
        for t in range(t_steps):
            for z in range(n_spatial):
                MassFluxDev[z, t] = abs(mdot_0 - (self.u_res[z, t] * self.reactor.rho_fl(self.w_i_res[z, t, :].T, self.T_res[z, t],
                                                                                    self.p_res[z, t]).__float__())) / mdot_0 * 100
        print("Maximal mass flux deviation: ", np.max(MassFluxDev))
