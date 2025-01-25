from V2.classes.Parameters.Component import Component
from V2.classes.Parameters.Discretization import Discretization

import casadi as CasADi
import numpy as np
import matplotlib.pyplot as plt


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

        # TODO make it work for 2D !
        # setting initial values for w_i
        for comp in range(n_components):
            if comp == 0:
                x_0[0: n_axial] = w_i_in[0]
            else:
                x_0[n_axial*comp : n_axial*(comp+1)] = w_i_in[comp]
        # setting initial values for T
        x_0[n_components * n_axial:] = T_in
        self.x_0 = x_0

        z_0= np.zeros(n_axial * n_radial * 2)
        z_0[0:n_axial * n_radial] = u_in
        z_0[n_axial * n_radial:] = p_in
        self.z_0 = z_0


########################################################################################################################

    def integrate(self):
        self.results =  self.integrator(x0=self.x_0, z0=self.z_0)

        #TODO below here only for first testing!

        n_axial = self.axial_discretization.num_volumes
        n_radial = self.radial_discretization.num_volumes
        n_components = self.reactor.n_components
        t_steps = self.time_discretization.num_volumes+1

        res_x = self.results['xf'].full()

        w_i_res = np.empty(shape=(n_axial, t_steps, n_components))
        T_res = np.empty(shape=(n_axial, t_steps))
        for t in range(t_steps):
            T_res[:, t] = res_x[n_components * n_axial:, t]
            for comp in range(n_components):
                w_i_res[:, t, comp] = res_x[comp * n_axial: n_axial * (comp + 1), t]

        ae_res = self.results['zf'].full()

        u_res = ae_res[:n_axial, :]
        p_res = ae_res[n_axial:, :]

        w_i_in = self.reactor.w_i_in
        T_in = self.reactor.T_in
        u_in = self.reactor.u_in
        p_in = self.reactor.p_in


        MassFluxDev = np.empty(shape=(n_axial, t_steps))
        mdot_0 = (u_in * self.reactor.rho_fl(w_i_in, T_in, p_in)).__float__()
        for t in range(t_steps):
            for z in range(n_axial):
                MassFluxDev[z, t] = abs(mdot_0 - (u_res[z, t] * self.reactor.rho_fl(w_i_res[z, t, :].T, T_res[z, t],
                                                                      p_res[z, t]).__float__())) / mdot_0 * 100
                #print(MassFluxDev[z, t])
        print("Maximal mass flux deviation: ", np.max(MassFluxDev))

        w_f = CasADi.SX.sym('w_f', 4)
        T_f = CasADi.SX.sym('T_f')
        p_f = CasADi.SX.sym('p_f')
        u_f = CasADi.SX.sym('u_f')

        # ETA
        f_eta_casADI = CasADi.Function('f_eta_casADi', [w_f, T_f, p_f], [self.reactor.effFactor(w_f, T_f, p_f)])
        eta = np.empty(shape=(n_axial, t_steps))

        # check function
        f_check_casADI = CasADi.Function('f_check_casADI', [w_f, T_f, p_f, u_f], [])
        for t in range(t_steps):
            for z in range(n_axial):
                if z == 0:
                    pass
                eta[z, t] = f_eta_casADI(w_i_res[z, t, :], T_res[z, t], p_res[z, t])
                #print(eta[z, t])
                check = f_check_casADI(w_i_res[z, t, :], T_res[z, t], p_res[z, t], u_res[z, t])
                #print("w_i", w_i_res[z, t, :])
                #print("conv heat flux", self.reactor.axialMassFlow(T_res[z, t], w_i_res[z, t, :], w_i_res[z-1, t, :], u_res[z, t], p_res[z, t], 2))
                #print("rate eq.", check)



        # Plot results
        # Generate plot
        fig, axs = plt.subplots(5, 1, figsize=(4.2, 6.7), constrained_layout=True, sharex=True)
        # Define colors
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))
        # Plot
        t_step = 10

        axs[0].plot(w_i_res[:, t_step, 0], color=colors[0])
        axs[0].plot(w_i_res[:, t_step, 1], color=colors[1])
        axs[0].plot(w_i_res[:, t_step, 2], color=colors[2])
        axs[0].plot(w_i_res[:, t_step, 3], color=colors[3])
        axs[1].plot(T_res[:, t_step], color=colors[0])
        axs[1].set_ylim([300, 350])
        axs[2].plot(u_res[:, t_step], color=colors[0])
        axs[3].plot(p_res[:, t_step], color=colors[0])

        axs[4].plot(eta[:, t_step], color=colors[0])
        axs[4].set_ylim([0, 1.1])
        # Axis
        axs[0].set_ylabel(r'$w_{\mathregular{A}}$')
        axs[1].set_ylabel(r'$T \; \mathregular{/K}$')
        axs[2].set_ylabel(r'$u \; \mathregular{/ms^{-1}}$')
        axs[3].set_ylabel(r'$p / Pa$')
        axs[4].set_ylabel(r'$eff Factor / -$')
        axs[4].set_xlabel(r'$z/L$')
        plt.show()