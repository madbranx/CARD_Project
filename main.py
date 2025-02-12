
from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator
import numpy as np

# DISCRETIZATION SETTINGS #TODO (TBD: DISCRETIZATION STUDY)

n_axial = 40
n_radial = 6
time_end = 3000
time_steps = 600
precision = 1e-10


## 1D CASE WITH VALIDATION PLOT

# reactor = FixedBedReactor(1, n_axial)
# reactor.setup()
# integrator = Integrator(reactor)
# integrator.setup(precision, precision, 0, time_end, time_steps)
# results = integrator.integrate()
#postprocessor = Postprocessor(reactor, "../results/01")
#postprocessor.plot_1D_vs_ValidationData("test", results, time_end)


## 2D CASE WITH RESULT PLOTS FOR T, P, U AND CONVERSION X_CO2
## WITH VALIDATION OF WALL TEMPERATURE
## WITH COMPARISION OF CENTER VALUES 1D VS 2D #TODO (TBD!)

reactor = FixedBedReactor(2, n_axial, n_radial)
reactor.setup()
integrator = Integrator(reactor)
integrator.setup(precision, precision, 0, time_end, time_steps)
results = integrator.integrate()
postprocessor = Postprocessor(reactor, "../results/02")
postprocessor.plot2D_Temperature("test2", results,time_steps)
postprocessor.plot_Twall_vs_validation("test2", results, time_steps)


## 2D CASE WITH EXTINCTION AND IGNITION ARCS PLOTS

# # setting up reactor and integrator
# reactor = FixedBedReactor(2, n_axial, n_radial)
# reactor.T_wall = 550 # Temperature for ignited reactor
# reactor.setup()
# integrator = Integrator(reactor)
# integrator.setup(precision, precision, 0, time_end, time_steps)
#
# # get results of ignited reactor
# result_ignited = integrator.integrate()
# w_i, T, p, u = result_ignited.get_rawValues()
#
# # calculating arcs
# results_ignition = []
# results_extinction = []
# T_walls = np.linspace(400, 550, 25)
# for T_wall in T_walls:
#     print("T_wall = ", T_wall)
#     # calculating ignition arcs
#     reactor.T_wall = T_wall
#     reactor.setup()
#     integrator.refresh()
#     result_ignition = integrator.integrate()
#     results_ignition.append(result_ignition)
#     # calculating extinction arcs
#     integrator.refresh()
#     integrator.set_specific_InitialValues(w_i, T, p, u)
#     result_extinction = integrator.integrate()
#     results_extinction.append(result_extinction)
#
# postprocessor = Postprocessor(reactor, "../results/02")
# postprocessor.plot_ignitionArc(results_ignition, results_extinction, T_walls, time_steps)


# # TESTING FUNCTION
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


# Create Ignition/Extinction Arcs:
#   calc till steady state, plot X over T wall
#   extinction, last steady state x=1 as starting values


#TODO
# Fix numerical issues & radial dispersion
# .
# .
# Create ULM diagram             TBD
# add Diskretisierungsformeln in PDF
# .
# Postprocessing (2D)            WORK IN PROGRESS
# .     Extinction / Ignition ARC Plots
# .     Effectiveness etc. Plot
# .     RMSE/MAPE for one Value
# .     Center 2D vs 1D calc / vs 1D validation data
# .
