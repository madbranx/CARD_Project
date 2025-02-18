from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator
import numpy as np


class Studies:
    def __init__(self):
        pass

    '''####################################### Validation #######################################'''

    def validation(self, foldername, n_axial, n_radial, time_end, timesteps_1D, timesteps_2D, log=False):

        result_1D = self.__validation_1D(n_axial, time_end, timesteps_1D, log=log)
        result_2D = self.__validation_2D(n_axial, n_radial, time_end, timesteps_2D, log=log)

        postprocessor = Postprocessor("results")
        postprocessor.plot_1D_and_2D_vs_ValidationData(foldername, result_1D, result_2D, timesteps_1D, timesteps_2D)


    def __validation_1D(self, n_axial, time_end, time_steps, t_equi=False, z_equi=False, log=False):
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

        return results


    def __validation_2D(self, n_axial, n_radial, time_end, time_steps, t_equi=False, z_equi=False, r_equi=False, log=False):
        ## 2D CASE WITH VALIDATION PLOT
        options_int = {
            "tols": [1e-4, 1e-4],
            "get_runtime": True,
            "log": log,
        }

        reactor = FixedBedReactor(2, n_axial, n_radial, z_equi=z_equi, r_equi=r_equi)
        reactor.setup()

        integrator = Integrator(reactor)
        integrator.set_options(**options_int)
        integrator.setup(0, time_end, time_steps, t_equi=t_equi)
        results = integrator.integrate()

        return results

    '''######################################## Base Case #######################################'''

    def base_case_2D(self, foldername, n_axial, n_radial, time_end, time_steps, t_equi=False, z_equi=False, r_equi=False, log=False):
        options_int = {
            "tols": [1e-4, 1e-4],
            "get_runtime": True,
            "log": log,
        }

        reactor = FixedBedReactor(2, n_axial, n_radial, z_equi=z_equi, r_equi=r_equi)
        reactor.setup()

        integrator = Integrator(reactor)
        integrator.set_options(**options_int)
        integrator.setup(0, time_end, time_steps, t_equi=t_equi)
        results = integrator.integrate()

        postprocessor = Postprocessor("results")
        postprocessor.plot2D_wTpu_X(foldername, results, time_steps, 2) # 2 = CO2

    '''##################################### Discretization #####################################'''

    def discretization_study(self, foldername, time_end, time_steps_axial, n_axial_ref, n_axials, time_steps_radial, n_axial_radial, n_radial_ref, n_radials, log=False):
        # Setting Simulation Options
        options_axial = {
            "tols": [1e-8, 1e-8],
            "get_runtime": False,
            "log": log,
        }

        options_radial = {
            "tols": [1e-4, 1e-4],
            "get_runtime": False,
            "log": log,
        }

        def run_sim(reactor, options, timesteps):
            reactor.setup()
            integrator = Integrator(reactor)
            integrator.set_options(**options)
            integrator.setup(0, time_end, timesteps)
            return integrator.integrate()

        def print_progress(dim, n_ax, n_rad):
            print("simulating " + str(dim) + "D with # ax|rad = " + str(n_ax) + "|" + str(n_rad))

        ## 1) Axial Simulations
        # Simulating axial reference
        reactor = FixedBedReactor(1, n_axial_ref, 1)
        print_progress(1, n_axial_ref, 1)
        result_ax_ref = run_sim(reactor, options_axial, time_steps_axial)

        # Simulating axial ed & ned
        results_ax_ed = []
        results_ax_ned = []
        for n_axial in n_axials:
            print_progress(1, n_axial, 1)
            # ed
            reactor = FixedBedReactor(1, n_axial, 1, z_equi=True)
            results_ax_ed.append(run_sim(reactor, options_axial, time_steps_axial))
            # ned
            reactor = FixedBedReactor(1, n_axial, 1)
            results_ax_ned.append(run_sim(reactor,options_axial,time_steps_axial))

        ## 1) Radial Simulations
        # Simulating raidal reference
        reactor = FixedBedReactor(2, n_axial_radial, n_radial_ref)
        print_progress(2, n_axial_radial, n_radial_ref)
        result_rad_ref = run_sim(reactor, options_radial, time_steps_radial)

        # Simulating radial ed & ned
        results_rad_ed = []
        results_rad_ned = []
        for n_radial in n_radials:
            print_progress(2, n_axial_radial, n_radial)
            # ed
            reactor = FixedBedReactor(2, n_axial_radial, n_radial, z_equi=True, r_equi=True)
            results_rad_ed.append(run_sim(reactor, options_radial, time_steps_radial))
            # ned
            reactor = FixedBedReactor(2, n_axial_radial, n_radial)
            results_rad_ned.append(run_sim(reactor, options_radial, time_steps_radial))

        ## Postprocessing
        postprocessor = Postprocessor("results")
        postprocessor.plot_disdiscretizationStudy(foldername, result_ax_ref, results_ax_ed, results_ax_ned, result_rad_ref, results_rad_ed, results_rad_ned, time_steps_axial, time_steps_radial)






# TODO: update methods below with new syntax!





        #
        # plot_times = [1, 5, 10, 50, 100]  # in %
        # for plot_time in plot_times:
        #     postprocessor.plot2D_wTpu_X("test2", results, int(plot_time / 100 * time_steps))
        #




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

