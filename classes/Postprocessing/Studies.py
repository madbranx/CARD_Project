from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator
import numpy as np
import time

"""
In the Studies class, all the realized studies of the reactor are defined in methods. 
If a studie is to be carried out, respective methode has to be called with the required parameters in the main.py file.
It will automatically set up the FixedBedReactor and Integrator class and call them with the settings of the studie.
When the studie is finished and the results acquired, the respective postprocessor methode is called and the results evaluated.
"""

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

    '''######################################## Base Case ########################################'''

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

    '''##################################### Discretization ######################################'''

    def discretization_study(self, foldername, time_end, time_steps_axial, n_axial_ref, n_axials, time_steps_radial, n_axial_radial_ref, n_radial_ref, n_radials, log=False):
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
            'max_step_size': 0.1,
            "max_num_steps": 20000,
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
        # Simulating radial reference
        reactor = FixedBedReactor(2, n_axial_radial_ref, n_radial_ref)
        print_progress(2, n_axial_radial_ref, n_radial_ref)
        result_rad_ref = run_sim(reactor, options_radial, time_steps_radial)

        # Simulating radial ed & ned
        results_rad_ed = []
        results_rad_ned = []
        for n_radial in n_radials:
            n_axial = n_axial_radial_ref # = int(np.ceil(n_radial*n_axial_radial_ref/n_radial_ref))
            print_progress(2, n_axial, n_radial)
            # ed
            reactor = FixedBedReactor(2, n_axial, n_radial, z_equi=True, r_equi=True)
            results_rad_ed.append(run_sim(reactor, options_radial, time_steps_radial))
            # ned
            reactor = FixedBedReactor(2, n_axial, n_radial)
            results_rad_ned.append(run_sim(reactor, options_radial, time_steps_radial))

        ## Postprocessing
        postprocessor = Postprocessor("results")
        postprocessor.plot_disdiscretizationStudy(foldername, result_ax_ref, results_ax_ed, results_ax_ned, result_rad_ref, results_rad_ed, results_rad_ned, time_steps_axial, time_steps_radial)

    def measure_time_discretization(self, n_axial, n_radial, t_steps, t_end):
        import os
        os.environ['CASADI_LOG_LEVEL'] = 'ERROR'

        print("#ax = ", n_axial)
        print("#rad = ", n_radial)
        print("#t-steps = ", t_steps)
        print("sim time = ", t_end, "\n\n")

        options = {
            "tols": [1e-4, 1e-4],
            "get_runtime": False,
            'max_step_size': 1,
            "max_num_steps": 20000,
        }

        def time_sim(z_eqi, r_equi, t_equi):
            reactor = FixedBedReactor(2, n_axial, n_radial, z_equi=z_eqi, r_equi=r_equi)
            reactor.setup()
            integrator = Integrator(reactor)
            integrator.set_options(**options)
            integrator.setup(0, t_end, t_steps, t_equi=t_equi)

            start_time = time.time()  # Start timer

            integrator.integrate()

            end_time = time.time()  # End timer

            z, r, t = "x", "x", "x"
            if z_eqi: z = "o"
            if r_equi: r = "o"
            if t_equi: t = "o"
            print(z + "|" + r + "|" + t + f"    {end_time - start_time:.1f} seconds")

        print("x = non equi-distand, o = equi-distant\nz|r|t    time      o = non equi-distand, x = equi-distant\n")
        time_sim(True, True, True)
        time_sim(True, False, True)
        time_sim(False, True, True)
        time_sim(False, False, True)
        time_sim(True, True, False)
        time_sim(False, False, False)

    '''################################# Ignition/Extinction Arcs ################################'''

    def arcs(self, foldername, dim, time_end, time_steps, T_walls,n_axial, n_radial=1, log=False, plotting = True):
        ## 2D EXTINCTION AND IGNITION ARCS PLOTS

        # Setting Simulation Options
        if dim == 1:
            n_radial = 1
            options = {
                "tols": [1e-8, 1e-8],
                "get_runtime": False,
                "log": log,
            }
        else:
            options = {
                "tols": [1e-4, 1e-4],
                "get_runtime": False,
                "log": log,
                'max_step_size': 0.1,
                "max_num_steps": 20000,
            }

        def run_sim(reactor, w_i=None, T=None, p=None, u=None):
            reactor.setup()
            integrator = Integrator(reactor)
            integrator.set_options(**options)
            integrator.setup(0, time_end, time_steps)
            if w_i is not None:
                integrator.set_specific_InitialValues(w_i, T, p, u)
            return integrator.integrate()

        # Setting up Reactor
        reactor = FixedBedReactor(dim, n_axial, n_radial, z_equi=True)
        w_i, T, p, u = None, None, None, None

        # # Simulating extinguished reactor -> not needed currently
        # reactor.T_wall = T_walls[0]
        # result_extinguished = run_sim(reactor)
        # w_i, T, p, u = result_extinguished.get_rawValues(time_steps)

        # calculating arcs
        results_ignition = []
        for T_wall in T_walls:

            # setting T_wall
            print("Ignition: T_wall = ", T_wall)
            reactor.T_wall = T_wall
            reactor.setup()
            try: # running simulation and setting new starting values
                results_ignition.append(run_sim(reactor, w_i, T, p, u))
                w_i, T, p, u  = results_ignition[-1].get_rawValues(time_steps)
            except:
                print("ignition failed")
                pass

        # extinction arcs
        T_walls_ext = np.flip(T_walls)
        results_extinction = []
        for T_wall in T_walls_ext:

            # setting T_wall
            print("Extinction: T_wall = ", T_wall)
            reactor.T_wall = T_wall
            reactor.setup()
            try: # running simulation and setting new starting values
                results_extinction.append(run_sim(reactor, w_i, T, p, u))
                w_i, T, p, u = results_extinction[-1].get_rawValues(time_steps)
            except:
                print("extinction failed")
                pass

        if plotting == False:
            return results_ignition, results_extinction, T_walls, T_walls_ext
        else:
            ## Postprocessing
            postprocessor = Postprocessor("results")
            postprocessor.plot_ignitionArc(foldername, results_ignition, results_extinction, T_walls, T_walls_ext, time_steps)

    def combinedArcs(self, foldername, time_end, time_steps, T_walls, n_axial, n_radial, log=False):
        results_ignition_1d, results_extinction_1d, T_walls_1d, T_walls_ext_1d = self.arcs(foldername, 1, time_end,
                                                                                           time_steps, T_walls,
                                                                                           n_axial * 2, log=log,
                                                                                           plotting=False)
        results_ignition_2d, results_extinction_2d, T_walls_2d, T_walls_ext_2d = self.arcs(foldername, 2, time_end,
                                                                                           time_steps, T_walls,
                                                                                           n_axial, n_radial,
                                                                                           log=log, plotting=False)

        ## Postprocessing
        postprocessor = Postprocessor("results")
        postprocessor.plot_ignitionArc_1_2D(foldername, results_ignition_1d, results_ignition_2d,
                                            results_extinction_1d, results_extinction_2d, T_walls_1d,
                                            T_walls_ext_1d, T_walls_2d, T_walls_ext_2d, time_steps)

        '''#################################### Catalyst Variation ###################################'''

    def cat_variation_diameter(self, foldername, n_axial, n_radial, time_end, time_steps, d_cats, log=False):
        options = {
            "tols": [1e-4, 1e-4],
            "get_runtime": False,
            "log": log,
            'max_step_size': 0.1,
            "max_num_steps": 20000,
        }

        postprocessor = Postprocessor("results")

        def run_sim(reactor):
            reactor.setup()
            integrator = Integrator(reactor)
            integrator.set_options(**options)
            integrator.setup(0, time_end, time_steps)
            return integrator.integrate()

        # Setting up Reactor
        reactor = FixedBedReactor(2, n_axial, n_radial, z_equi=True)

        # getting reference sim result for pressure drop
        result_ref = run_sim(reactor)


        # calculating d_cat variation
        results_d_cat = []
        for d_cat in d_cats:
            # setting d_cat
            print("cat diameter = ", d_cat)
            reactor.cat_diameter = d_cat
            reactor.eps = reactor.calculate_void_fraction()
            print("eps = ", reactor.eps)
            reactor.setup()
            try: # running simulation and setting new starting values
                results_d_cat.append(run_sim(reactor))
            except:
                print("simulation failed")
                pass

        postprocessor.plotCatVariation_diameter(foldername, "cat_diameter", result_ref, results_d_cat, d_cats * 1e3, time_steps)

    def cat_variation_pore(self, foldername, n_axial, n_radial, time_end, time_steps, d_pores, log=False):
        options = {
            "tols": [1e-4, 1e-4],
            "get_runtime": False,
            "log": log,
            'max_step_size': 0.1,
            "max_num_steps": 20000,
        }

        postprocessor = Postprocessor("results")

        def run_sim(reactor):
            reactor.setup()
            integrator = Integrator(reactor)
            integrator.set_options(**options)
            integrator.setup(0, time_end, time_steps)
            return integrator.integrate()

        # Setting up Reactor
        reactor = FixedBedReactor(2, n_axial, n_radial, z_equi=True)

        # getting reference sim result for pressure drop
        result_ref = run_sim(reactor)

        # calculating d_pores variation
        results_d_pores = []
        for d_pore in d_pores:
            # setting d_pore
            print("pore diameter = ", d_pore)
            reactor.diameter_pore = d_pore
            reactor.setup()
            try: # running simulation and setting new starting values
                results_d_pores.append(run_sim(reactor))
            except:
                print("simulation failed")
                pass

        postprocessor.plotCatVariation_pore(foldername, "pore_diameter", result_ref, results_d_pores, d_pores * 1e9, time_steps)

    '''#################################### Ignition/Extinction behavior ###################################'''

    def ignition_behavior(self, foldername, n_axial, n_radial, time_end, timesteps, t_evals, threshold_T, threshold_X, log):
        options = {
            "tols": [1e-4, 1e-4],
            "get_runtime": False,
            "log": log,
            'max_step_size': 0.1,
            "max_num_steps": 10000,
        }
        reactor = FixedBedReactor(2, n_axial, n_radial)
        reactor.setup()
        integrator = Integrator(reactor)
        integrator.set_options(**options)
        integrator.setup(0, time_end, timesteps)
        result =  integrator.integrate()

        postprocessor = Postprocessor("results")
        postprocessor.plot_ignition_behavior(foldername, result, timesteps, time_end, t_evals, threshold_T,threshold_X)
