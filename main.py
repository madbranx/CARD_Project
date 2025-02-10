from classes.Postprocessing.Postprocessor import Postprocessor
from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator

reactor = FixedBedReactor(2, 30, 15)

integrator = Integrator(reactor)
integrator.setup(1e-5, 1e-5, 0, 1000, 100)
results = integrator.integrate()

#postprocessor = Postprocessor(reactor, "../results/01")
#postprocessor.plot_1D_vs_ValidationData("test", results, 100)

postprocessor = Postprocessor(reactor, "../results/02")
postprocessor.plot2D_Temperature("test2", results,10)
postprocessor.plot2D_Temperature("test2", results,100)


# Questions
# Create Ignition/Extinction Arcs after 2D implementation? - 1D and 2D
#   calc Till Kasselmann steady state, plot X over T wall
#   extinction, last steady state x=1 as starting values


#TODO V2
# add collision area for H2O
# .
# .
# Create ULM diagram             TBD
# add Log?
# .
# Conservations: Implement physics as CasADi functions (2D)
#       Species Conservation     WORK IN PROGRESS
#       Energy Conservation      WORK IN PROGRESS
# .
# Postprocessing (2D)            WORK IN PROGRESS
# .     2D Plots -> plot for constant value curves (?)
# .     Extinction / Ignition ARC Plots
