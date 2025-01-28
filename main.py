from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator

reactor = FixedBedReactor(1, 100)

integrator = Integrator(reactor)
integrator.setup(1e-11, 1e-11, 0, 1000, 100)
results = integrator.integrate()

postprocessor = Postprocessor(reactor, "../results/01")
postprocessor.plot_1D_vs_ValidationData("test", results, 100)


# import casadi as CasADi
# from classes.Parameters.Component import Component
# w = CasADi.SX([0, 0, 0.8, 0.2])
#
# print(reactor.massFraction_weighted_average(w, Component.HEAT_CAPACITY, 500))

# Questions
# Validation Data @ which timepoints/physics/material properties etc.?
# What causes difference to validation data? reaction rate?
# Create Ignition/Extinction Arcs after 2D implementation?
# Radial mass flow: how to get the second derivative?


#TODO V2
# add collision area for H2O
# ADD INPUT / OUTPUT UNITS FOR ALL FUNCTIONS as comments
# .
# Postprocessing                WORK IN PROGRESS
# .     Add Plotting variable Reactor method as Casadi function
# .     Add 1D plot without Validation, add 1D plot with 2 sim results, add difference plot for 2 sim
# .     Extinction / Ignition ARC Plots
# .
# Validate 1D Model             WORK IN PROGRESS
# .     Check Reaction Rate and Eff Factor etc.
# .
# .
# Create ULM diagram           TBD
# add Log?
# .
# Conservations: Implement physics as CasADi functions (2D)
#       Species Conservation     TBD
#       Energy Conservation      TBD
# Postprocessing (2D)
