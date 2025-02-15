from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator
import numpy as np

### SCROLL DOWN FOR MAIN FUNCTION CALLS ###

## SIMULATION FUNCTIONS:
###########################################################################################
def validation_1D(n_axial, time_end, time_steps):
    ## 1D CASE WITH VALIDATION PLOT

    reactor = FixedBedReactor(1, n_axial)
    reactor.T_wall = 550
    reactor.setup()

    integrator = Integrator(reactor)
    integrator.setup(0, time_end, time_steps)
    results = integrator.integrate()

    postprocessor = Postprocessor(reactor, "../results/1D_val")
    postprocessor.plot_1D_vs_ValidationData("test", results, time_steps)


def plot_2D_with_TwallVal(n_axial, n_radial, time_end, time_steps):
    ## 2D CASE WITH RESULT PLOTS FOR T, P, U AND CONVERSION X_CO2
    ## WITH VALIDATION OF TEMPERATURE at outermost and innermost cell

    reactor = FixedBedReactor(2, n_axial, n_radial)
    reactor.setup()
    integrator = Integrator(reactor)
    integrator.setup(0, time_end, time_steps)
    results = integrator.integrate()

    postprocessor = Postprocessor(reactor, "../results/2D_withTwallVal")
    plot_times = [1, 5, 10, 50, 100] # in %
    for plot_time in plot_times:
        postprocessor.plot2D_Temperature("test2", results, int(plot_time/100*time_steps))

    postprocessor.plot_Twall_vs_validation("test2", results, time_steps)


def pseudo_2D_vs_1D(n_axial, n_radial, time_end, time_steps):
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

def arcs_1d(n_axial, time_end, time_steps, ignited_T, T_walls):
    ## 1D CASE WITH EXTINCTION AND IGNITION ARCS PLOTS

    #setting up reactor and integrator
    reactor = FixedBedReactor(1, n_axial)
    reactor.T_wall = ignited_T # Temperature for ignited reactor
    reactor.setup()
    integrator = Integrator(reactor)

    # get results of ignited reactor
    integrator.setup(0, time_end, time_steps)
    result_ignited = integrator.integrate()
    postprocessor = Postprocessor(reactor, "../results/arcs_1d")
    w_i, T, p, u = result_ignited.get_rawValues()

    # calculating arcs
    results_ignition = []
    results_extinction = []
    for T_wall in T_walls:

        # calculating ignition arcs
        print("T_wall = ", T_wall)
        reactor.T_wall = T_wall
        reactor.setup()
        try:
            integrator.setup(0, time_end, time_steps)
            result_ignition = integrator.integrate()
            results_ignition.append(result_ignition)
        except:
            pass

        #extinction arcs
        try:
            integrator.set_specific_InitialValues(w_i, T, p, u)
            result_extinction = integrator.integrate()
            results_extinction.append(result_extinction)
        except:
            pass

    postprocessor.plot_ignitionArc1D(results_ignition, results_extinction, T_walls, time_steps)


def arcs_2d(n_axial, n_radial, time_end, time_steps, ignited_T, T_walls, plotting=False):
    ## 2D CASE WITH EXTINCTION AND IGNITION ARCS PLOTS

    # setting up reactor and integrator
    reactor = FixedBedReactor(2, n_axial, n_radial)
    reactor.T_wall = ignited_T # Temperature for ignited reactor
    reactor.setup()
    integrator = Integrator(reactor)

    # get results of ignited reactor
    integrator.setup(0, time_end, time_steps)
    result_ignited = integrator.integrate()
    postprocessor = Postprocessor(reactor, "../results/arcs_2d")
    w_i, T, p, u = result_ignited.get_rawValues()

    # calculating arcs
    results_ignition = []
    results_extinction = []
    done_T_walls = []
    for T_wall in T_walls:

        # calculating ignition arcs
        print("T_wall = ", T_wall)
        reactor.T_wall = T_wall
        reactor.setup()
        try:
            integrator.setup(0, time_end, time_steps)
            result_ignition = integrator.integrate()
        except:
            continue

        #extinction arcs
        try:
            integrator.setup(0, time_end, time_steps)
            integrator.set_specific_InitialValues(w_i, T, p, u)
            result_extinction = integrator.integrate()
            results_extinction.append(result_extinction)

            results_ignition.append(result_ignition)
            done_T_walls.append(T_wall)
            if plotting: # Plot each succesfull T_wall step
                postprocessor.plot_ignitionArc(results_ignition, results_extinction, done_T_walls, time_steps)
        except:
            pass

    # plot all T_wall steps at end
    postprocessor.plot_ignitionArc(results_ignition, results_extinction, T_walls, time_steps)


###########################################################################################

#validation_1D(40, 3000, 200)

#arcs_1d(60, 4000, 500, 550, np.linspace(350, 500, 5))

#plot_2D_with_TwallVal(150, 12, 3000, 6000)

#pseudo_2D_vs_1D(20, 5, 3000, 300)

arcs_2d(150, 12, 3000, 12000, 550, np.linspace(300, 550, 25), True)

###########################################################################################

#TODO
# Fix numerical issues & radial dispersion
# . -> add Flux limiter?/if/else for differences?
# .
# Create ULM diagram             TBD
# add Diskretisierungsformeln in PDF
# .
# Postprocessing (2D)            WORK IN PROGRESS
#       (TBD: DISCRETIZATION STUDY)
# .     Extinction / Ignition ARC Plots
# .     Effectiveness etc. Plot
# .     RMSE/MAPE for one Value

###########################################################################################

## TESTING FUNCTIONS

# import casadi as CasADi
# T = 500
# p = 5e5
# u = 1
#
# w1 = CasADi.SX([0.25, 0.25, 0.25, 0.25])
# w2 = CasADi.SX([0, 0, 0.800001, 0.199999])
#
# wTpu1 = [w1, T, p, u]
# wTpu2 = [w2, T, p, u]
#
# print("j dispersion = ", reactor.calc_j_dispersion(0.1, 0.1, 0.1, wTpu1, wTpu2, 2))
#
# print("lambda mass avg = ",reactor.massFraction_weighted_average(w1, Component.THERMAL_CONDUCTIVITY, T))
# print("lambda mixture rule = ",reactor.calc_fluid_conductivity(T, w1))
#
# print("wall contact HT coeff = ", reactor.calc_heatTransferCoefficient_contact(T, p , w1))
# print("k reactor wall / R = ", reactor.calc_resistanceWall()/(reactor.reactorDiameter/2))
# print(reactor.calc_eff_disp_coeff(wTpu1, 3))



