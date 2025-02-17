from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator
import numpy as np


class Studies:
    def __init__(self):
        pass

    def validation_1D(self, n_axial, time_end, time_steps, t_equi=False, z_equi=False, log=False):
        ## 1D CASE WITH VALIDATION PLOT
        options_int = {
            "tols": [1e-10, 1e-10],
            "get_runtime": True,
            "log": log,
        }

        reactor = FixedBedReactor(1, n_axial, z_equi=z_equi)
        reactor.setup()

        integrator = Integrator(reactor)
        integrator.set_options(**options_int)
        integrator.setup(0, time_end, time_steps, t_equi=t_equi)
        results = integrator.integrate()

        postprocessor = Postprocessor(reactor, "../results/1D_val")
        postprocessor.plot_1D_vs_ValidationData("test", results, time_steps)

    # TODO: update methods below with new syntax!

    def plot_2D_with_TwallVal(self, n_axial, n_radial, time_end, time_steps):
        ## 2D CASE WITH RESULT PLOTS FOR T, P, U AND CONVERSION X_CO2
        ## WITH VALIDATION OF TEMPERATURE at outermost and innermost cell

        reactor = FixedBedReactor(2, n_axial, n_radial)
        reactor.setup()
        integrator = Integrator(reactor)
        integrator.setup(0, time_end, time_steps)
        results = integrator.integrate()

        postprocessor = Postprocessor(reactor, "../results/2D_withTwallVal")
        plot_times = [1, 5, 10, 50, 100]  # in %
        for plot_time in plot_times:
            postprocessor.plot2D_wTpu_X("test2", results, int(plot_time / 100 * time_steps))

        postprocessor.plot_Twall_vs_validation("test2", results, time_steps)

    def pseudo_2D_vs_1D(self, n_axial, n_radial, time_end, time_steps):
        # 2D with 1 axial element
        reactor_2D = FixedBedReactor(2, n_axial, n_radial)
        reactor_2D.setup()
        integrator = Integrator(reactor_2D)
        integrator.setup(0, time_end, time_steps)
        results_2D = integrator.integrate()

        # 1D
        reactor_1D = FixedBedReactor(1, n_axial)
        reactor_1D.setup()
        integrator = Integrator(reactor_1D)
        integrator.setup(0, time_end, time_steps)
        results_1D = integrator.integrate()

        # Plotting
        postprocessor = Postprocessor(reactor_2D, "../results/p2D_vs_1D")
        postprocessor.plot1D_vsPseudo1D_vs_val("test2", results_1D, results_2D, time_steps)

    def arcs_1d(self, n_axial, time_end, time_steps, T_walls):
        ## 1D CASE WITH EXTINCTION AND IGNITION ARCS PLOTS

        # setting up reactor and integrator
        reactor = FixedBedReactor(1, n_axial, 1)
        reactor.T_wall = T_walls[0]  # Temperature for extinguished reactor
        reactor.setup()
        integrator = Integrator(reactor)

        # get results of unignited reactor
        integrator.setup(0, time_end, time_steps)
        result_extinguished = integrator.integrate()
        postprocessor = Postprocessor(reactor, "../results/arcs_1d")
        w_i, T, p, u = result_extinguished.get_rawValues(time_steps)

        # calculating arcs
        results_ignition = []
        for T_wall in T_walls:

            # calculating ignition arcs
            print("Ignition: T_wall = ", T_wall)
            reactor.T_wall = T_wall
            reactor.setup()
            try:
                integrator.setup(0, time_end, time_steps)
                integrator.set_specific_InitialValues(w_i, T, p, u)
                result_ignition = integrator.integrate()
                results_ignition.append(result_ignition)
                w_i, T, p, u = result_ignition.get_rawValues(time_steps)
            except:
                pass

        # extinction arcs
        T_walls_ext = np.flip(T_walls)
        results_extinction = []
        for T_wall in T_walls_ext:
            print("Extinction: T_wall = ", T_wall)
            reactor.T_wall = T_wall
            reactor.setup()
            try:
                integrator.setup(0, time_end, time_steps)
                integrator.set_specific_InitialValues(w_i, T, p, u)
                result_extinction = integrator.integrate()
                results_extinction.append(result_extinction)
                w_i, T, p, u = result_extinction.get_rawValues()
            except:
                pass

        postprocessor.plot_ignitionArc1D(results_ignition, results_extinction, T_walls, time_steps, True)

    def arcs_2d(self, n_axial, n_radial, time_end, time_steps, T_walls):
        ## 2D CASE WITH EXTINCTION AND IGNITION ARCS PLOTS

        # setting up reactor and integrator
        reactor = FixedBedReactor(1, n_axial, n_radial)
        reactor.T_wall = T_walls[0]  # Temperature for extinguished reactor
        reactor.setup()
        integrator = Integrator(reactor)

        # get results of unignited reactor
        integrator.setup(0, time_end, time_steps)
        result_extinguished = integrator.integrate()
        postprocessor = Postprocessor(reactor, "../results/arcs_2d")
        w_i, T, p, u = result_extinguished.get_rawValues(time_steps)

        # calculating arcs
        results_ignition = []
        for T_wall in T_walls:

            # calculating ignition arcs
            print("Ignition: T_wall = ", T_wall)
            reactor.T_wall = T_wall
            reactor.setup()
            try:
                integrator.setup(0, time_end, time_steps)
                integrator.set_specific_InitialValues(w_i, T, p, u)
                result_ignition = integrator.integrate()
                results_ignition.append(result_ignition)
                w_i, T, p, u = result_ignition.get_rawValues(time_steps)
            except:
                pass

        # extinction arcs
        T_walls_ext = np.flip(T_walls)
        results_extinction = []
        for T_wall in T_walls_ext:
            print("Extinction: T_wall = ", T_wall)
            reactor.T_wall = T_wall
            reactor.setup()
            try:
                integrator.setup(0, time_end, time_steps)
                integrator.set_specific_InitialValues(w_i, T, p, u)
                result_extinction = integrator.integrate()
                results_extinction.append(result_extinction)
                w_i, T, p, u = result_extinction.get_rawValues()
            except:
                pass

        postprocessor.plot_ignitionArc2D(results_ignition, results_extinction, T_walls, T_walls_ext, time_steps)

    def test_2D_extinction(self, n_axial, n_radial, time_end, time_steps, T_wall_ignited, T_wall_extinction):
        # setting up reactor and integrator
        reactor = FixedBedReactor(2, n_axial, n_radial)
        reactor.T_wall = T_wall_ignited  # Temperature for ignited reactor
        reactor.setup()
        integrator = Integrator(reactor)
        postprocessor = Postprocessor(reactor, "../results/arcs_2d")

        # get results of ignited reactor
        integrator.setup(0, time_end, time_steps)
        result_ignited = integrator.integrate()
        w_i, T, p, u = result_ignited.get_rawValues(time_steps)

        postprocessor.plot2D_wTpu_X("name", result_ignited, time_steps)

        # simulating extinction
        reactor.T_wall = T_wall_extinction
        reactor.setup()
        integrator.setup(0, time_end, time_steps)
        integrator.set_specific_InitialValues(w_i, T, p, u)
        result_extinction = integrator.integrate()

        postprocessor.plot2D_wTpu_X("name", result_extinction, time_steps)

    def discretization_study1D(self, time_end, time_steps, n_axial_ref, n_axials):

        reactor = FixedBedReactor(1, n_axial_ref, 1)
        reactor.setup()
        integrator = Integrator(reactor)
        integrator.setup(0, time_end, time_steps)
        result_ref = integrator.integrate()

        results = []
        for n_axial in n_axials:
            reactor = FixedBedReactor(1, n_axial, 1)
            reactor.setup()
            integrator = Integrator(reactor)
            integrator.setup(0, time_end, time_steps)
            results.append(integrator.integrate())

        postprocessor = Postprocessor(reactor, "../results/discr_1d")
        postprocessor.plot_discretizationStudy1D(results, result_ref, time_steps)

    def discretization_study_radial(self, time_end, timesteps, n_radial_ref, n_axial, n_radials):

        reactor = FixedBedReactor(2, n_axial, n_radial_ref)
        reactor.setup()
        integrator = Integrator(reactor)
        integrator.setup(0, time_end, timesteps)
        result_ref = integrator.integrate()

        results = []
        for n_radial in n_radials:
            reactor = FixedBedReactor(2, n_axial, n_radial)
            reactor.setup()
            integrator = Integrator(reactor)
            integrator.setup(0, time_end, timesteps)
            results.append(integrator.integrate())

        postprocessor = Postprocessor(reactor, "../results/discr_radial")
        postprocessor.plot_discretizationStudy_radial(results, result_ref, timesteps)
