from classes.Discretization.Discretization import Discretization

import casadi as CasADi
import numpy as np
import matplotlib.pyplot as plt


class Integrator:

    def __init__(self, log, reactor):
        self.log = log
        log.addEntry("initializing Integrator", 0)
        self.reactor = reactor
        self.discretization = None
        self.integrator = None
        self.x_0 = None
        self.z_0 = None
        self.results = None

    def setup(self, abstol, reltol, t_start, t_stop, t_steps):
        self.log.addEntry("setting up Integrator", 0)
        self.log.addEntry("creating time discretization", 1)

        self.discretization = Discretization(self.log, t_steps, Discretization.EQUIDISTANT, start=t_start, end=t_stop)

        self.log.addEntry("tolerances", 1)
        self.log.addEntry("abstol = " + str(abstol), 2)
        self.log.addEntry("reltol = " + str(reltol), 2)

        options = {'abstol': abstol, 'reltol': reltol}
        dae = self.reactor.getDAEstruct()
        timepoints = self.discretization.get_faces()
        self.integrator = CasADi.integrator('I', 'idas', dae, timepoints[0], timepoints, options)

        self.__setInitialValues()

    def integrate(self):
        self.log.addEntry("staring to integrate", 0)
        self.results =  self.integrator(x0=self.x_0, z0=self.z_0)
        self.log.addEntry("finished to integrate", 0)

        #TODO below here only for first testing!

        n_axial, n_radial = self.reactor.getSpatialDiscretizations()
        n_comps = self.reactor.getNComponents()
        t_steps = self.discretization.num_faces

        res_x = self.results['xf'].full()

        w_i_res = np.empty(shape=(n_axial, t_steps, n_comps))
        T_res = np.empty(shape=(n_axial, t_steps))
        for t in range(t_steps):
            T_res[:, t] = res_x[n_comps * n_axial:, t]
            for comp in range(n_comps):
                w_i_res[:, t, comp] = res_x[comp * n_axial: n_axial * (comp + 1), t]

        ae_res = self.results['zf'].full()

        u_res = ae_res[:n_axial, :]
        p_res = ae_res[n_axial:, :]

        w_i_in, T_in, u_in, p_in = self.reactor.getInputValues()

        from classes.GeneralConversionFunctions.GeneralConversionFunctions import GeneralConversionFunctions
        RSQ = self.reactor.getRSQ()
        GCF = GeneralConversionFunctions(self.log, RSQ)

        MassFluxDev = np.empty(shape=(n_axial, t_steps))
        mdot_0 = (u_in * GCF.rho_fl(w_i_in, T_in, p_in)).__float__()
        for t in range(t_steps):
            for z in range(n_axial):
                MassFluxDev[z, t] = abs(mdot_0 - (u_res[z, t] * GCF.rho_fl(w_i_res[z, t, :].T, T_res[z, t],
                                                                      p_res[z, t]).__float__())) / mdot_0 * 100
                #print(MassFluxDev[z, t])
        print("Maximal mass flux deviation: ", np.max(MassFluxDev))

        w_f = CasADi.SX.sym('w_f', 4)
        T_f = CasADi.SX.sym('T_f')
        p_f = CasADi.SX.sym('p_f')
        u_f = CasADi.SX.sym('u_f')

        # ETA
        from classes.FixedBedReactor.SpeciesConservation.ChangeByReaction.EffectivenessFactor.EffectivenessFactor import \
            EffectivenessFactor
        eff_factor = EffectivenessFactor(self.log, RSQ, GCF)
        f_eta_casADI = CasADi.Function('f_eta_casADi', [w_f, T_f, p_f], [eff_factor.calc(w_f, T_f, p_f)])
        eta = np.empty(shape=(n_axial, t_steps))

        # AxMassFlow
        from classes.FixedBedReactor.SpeciesConservation.AxialMassFlow.AxialMassFlow import AxialMassFlow
        axMassFlow_inst = AxialMassFlow(self.log, GCF)
        #f_axMassFlow_CasADi = CasADi.Function('f_eta_casADi', [T_f,w_f, u_f, p_f], [axMassFlow_inst.calc(T_f,w_f, u_f, p_f, 2)])

        axMassFLow = np.empty(shape=(n_axial, t_steps))

        for t in range(t_steps):
            for z in range(n_axial):
                eta[z, t] = f_eta_casADI(w_i_res[z, t, :], T_res[z, t], p_res[z, t])
                #axMassFLow[z, t] = f_axMassFlow_CasADi( T_res[z, t], w_i_res[z, t, :],w_i_res[z-1, t, :], u_res[z, t], p_res[z, t])

                #print(axMassFLow[z, 100])


        # Plot results
        # Generate plot
        fig, axs = plt.subplots(5, 1, figsize=(4.2, 6.7), constrained_layout=True, sharex=True)
        # Define colors
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))
        # Plot
        t_step = 100

        axs[0].plot(w_i_res[:, t_step, 0], color=colors[0])
        axs[0].plot(w_i_res[:, t_step, 1], color=colors[1])
        axs[0].plot(w_i_res[:, t_step, 2], color=colors[2])
        axs[0].plot(w_i_res[:, t_step, 3], color=colors[3])
        axs[1].plot(T_res[:, t_step], color=colors[0])
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


    def __setInitialValues(self):
        self.log.addEntry("setting initial values", 1)
        w_i_in, T_in, u_in, p_in = self.reactor.getInputValues()

        n_axial, n_radial = self.reactor.getSpatialDiscretizations()
        n_comps = self.reactor.getNComponents()
        if n_radial is None:
            n_radial = 1

        x_0 = np.zeros(n_axial * n_radial * (n_comps + 1))

        # TODO make it work for 2D !
        # setting initial values for w_i
        for comp in range(n_comps):
            if comp == 0:
                x_0[0: n_axial] = w_i_in[0]
            else:
                x_0[n_axial*comp : n_axial*(comp+1)] = w_i_in[comp]
        # setting initial values for T
        x_0[n_comps * n_axial:] = T_in
        self.x_0 = x_0

        z_0= np.zeros(n_axial * n_radial * 2)
        z_0[0:n_axial * n_radial] = u_in
        z_0[n_axial * n_radial:] = p_in
        self.z_0 = z_0