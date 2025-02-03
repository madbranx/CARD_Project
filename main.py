from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator

reactor = FixedBedReactor(2, 50, 20)

integrator = Integrator(reactor)
integrator.setup(1e-5, 1e-5, 0, 1000, 100)
results = integrator.integrate()

#postprocessor = Postprocessor(reactor, "../results/01")
#postprocessor.plot_1D_vs_ValidationData("test", results, 100)

postprocessor = Postprocessor(reactor, "../results/02")
postprocessor.plot2D("test2", results,10)

# Questions
# Create Ignition/Extinction Arcs after 2D implementation? - 1D and 2D
#   calc Till Kasselmann steady state, plot X over T wall
#   extinction, last steady state x=1 as starting values


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
