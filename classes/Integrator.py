from classes.Parameters.Discretization import Discretization
from classes.Postprocessing.Results import Results

import casadi as CasADi
import numpy as np
import copy

"""
The Integrator class contains the setup of the solver as well as the solver settings and result extraction.
"""

class Integrator:
    def __init__(self, reactor):
        self.reactor = reactor

        self.time_discretization = None
        self.axial_discretization = reactor.axial_discretization
        self.radial_discretization = reactor.radial_discretization

        self.integrator = None
        self.options = None

        self.x_0 = None
        self.z_0 = None
        self.results = None

        self.w_i_res = None
        self.T_res = None
        self.u_res = None
        self.p_res = None

    def setup(self, t_start, t_stop, t_steps, **kwargs):
        if kwargs.get('t_equi', False) is True:
            self.time_discretization = Discretization(t_steps, start=t_start, end=t_stop)
        else:
            if kwargs.get('time_ranges', None) is not None: # when explicit array of time ranges are given
                ranges_t = kwargs['time_ranges']
            else: # otherwise determine ranges:
                if t_stop < 100: # equidistant for <100 s simulation time
                    ranges_t = [[t_start, t_stop, 1]]
                elif t_stop < 1000:
                    ranges_t = [[t_start, 90, 6], [90, t_stop, 3]]
                else: # non-equidistant with higher resolutions at the beginning of the simulation
                    ranges_t = [[t_start, 90, 6], [90, t_stop*0.1, 3], [t_stop*0.1, t_stop/3, 2], [t_stop/3, t_stop, 0.25]]

            self.time_discretization = Discretization(t_steps, Discretization.RELATIVE_ARRAY, ranges=ranges_t)

        self.refresh()


    def refresh(self):
        # Methode to change solver options.
        # Standard options are used if none were explicitly given before
        if self.options is None:
            self.set_options()

        dae = self.reactor.DAE
        timepoints = self.time_discretization.get_faces()

        self.integrator = CasADi.integrator('I', 'idas', dae, timepoints[0], timepoints, self.options)
        self.__setInitialValues()

    def set_options(self, **kwargs):
        # Methode with standart solver options. Use refresh methode to change settings from standart.
        [abstol, reltol] = kwargs.get('tols', [1e-8, 1e-8])
        if kwargs.get('log', False) is True:
            verbose = True
            disable_internal_warnings = False
        else:
            verbose = False
            disable_internal_warnings = True

        self.options = {
            'abstol': abstol,
            "scale_abstol": kwargs.get('scale_tol', False),
            'reltol': reltol,
            "step0": kwargs.get('init_step', 0.01),
            "max_step_size": kwargs.get('max_step_size', 1),
            "max_num_steps": kwargs.get('max_num_steps', 10000),
            "newton_scheme": kwargs.get('newton_scheme', 'direct'),
            "print_time": kwargs.get('get_runtime', False),
            "verbose": verbose,
            "disable_internal_warnings": disable_internal_warnings,
        }

    def set_specific_InitialValues(self, w_i, T, p, u ):
        # Methode to set specific initial conditions
        n_axial = self.axial_discretization.num_volumes
        n_radial = self.radial_discretization.num_volumes
        n_components = self.reactor.n_components

        x_0 = np.zeros(n_axial * n_radial * (n_components + 1))
        for comp in range(n_components):
            if comp == 0:
                x_0[0: n_axial * n_radial] = w_i[:, 0]
            else:
                x_0[n_axial * n_radial * comp: n_axial * n_radial * (comp + 1)] = w_i[:, comp]
        x_0[n_components * n_axial * n_radial:] = T[:]
        self.x_0 = x_0

        z_0 = np.zeros(n_axial * n_radial * 2)
        z_0[0:n_axial * n_radial] = u[:]
        z_0[n_axial * n_radial:] = p[:]
        self.z_0 = z_0

    def __setInitialValues(self):
        # Methode to set general initial conditions defined in the parameters class
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
        # Methode to start the integration with the set starting conditions and options
        self.results =  self.integrator(x0=self.x_0, z0=self.z_0)
        self.__extractResults()
        #self.__printMassDeviation()

        results = Results(self.axial_discretization, self.radial_discretization, self.time_discretization)
        results.set_reactor(self.reactor)
        results.add_values(copy.deepcopy(self.w_i_res), copy.deepcopy(self.T_res), copy.deepcopy(self.u_res), copy.deepcopy(self.p_res))
        return results

    def __extractResults(self):
        # Methode to convert results of integration into usable variables for the defined timesteps
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
        #print("Maximal mass flux deviation: ", np.max(MassFluxDev), "\n")
