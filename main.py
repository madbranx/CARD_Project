from classes.Parameters.Component import Component
from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator

reactor = FixedBedReactor(2, 30, 15)

integrator = Integrator(reactor)
integrator.setup(1e-8, 1e-8, 0, 5000, 250)
results = integrator.integrate()

#postprocessor = Postprocessor(reactor, "../results/01")
#postprocessor.plot_1D_vs_ValidationData("test", results, 100)

postprocessor = Postprocessor(reactor, "../results/02")
postprocessor.plot2D_Temperature("test2", results,10)
postprocessor.plot2D_Temperature("test2", results,100)
postprocessor.plot2D_Temperature("test2", results,250)

postprocessor.plot_Twall_vs_validation("test2", results, 100)

# # TESTING FUNCTION
import casadi as CasADi
T = 500
p = 5e5
u = 1

w1 = CasADi.SX([0.25, 0.25, 0.25, 0.25])
w2 = CasADi.SX([0, 0, 0.800001, 0.199999])

wTpu1 = [w1, T, p, u]
wTpu2 = [w2, T, p, u]

print("j dispersion = ", reactor.calc_j_dispersion(0.1, 0.1, 0.1, wTpu1, wTpu2, 2))

print("lambda mass avg = ",reactor.massFraction_weighted_average(w1, Component.THERMAL_CONDUCTIVITY, T))
print("lambda mixture rule = ",reactor.calc_fluid_conductivity(T, w1))

print("wall contact HT coeff = ", reactor.calc_heatTransferCoefficient_contact(T, p , w1))
print("k reactor wall / R = ", reactor.calc_resistanceWall()/(reactor.reactorDiameter/2))
print(reactor.calc_eff_disp_coeff(wTpu1, 3))


# Create Ignition/Extinction Arcs:
#   calc till steady state, plot X over T wall
#   extinction, last steady state x=1 as starting values


#TODO
# check alpha contact wall -> mean free path?
# fix radial mass flow
# add collision area for H2O
# lambda fl mixture vs. lambda mass frac
# .
# .
# Create ULM diagram             TBD
# Diskretisierungsformeln in PDF
# .
# Postprocessing (2D)            WORK IN PROGRESS
# .     Extinction / Ignition ARC Plots
# .     Effectiveness etc. Plot
# .     RMSE/MAPE for one Value
# .     Center 2D vs 1D calc / vs 1D validation data
# .
