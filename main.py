from pandas.plotting import plot_params

from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
from classes.Log.Log import Log

log = Log("first simulation")
reactor = FixedBedReactor(log, FixedBedReactor.ONE_D, 1)
reactor.setUp()
##
dae = reactor.getDAEstruct()
print(dae)

#TODO
# add all parameters and material properties to FixedBedReactor()
# Conservations: Implement physics as CasADi functions (1D)
#       Pressure drop            DONE
#       Mass conservation        DONE
#       Species Conservation     DONE
#       Energy Conservation      TBD
# .
# .
# Integrator
# Postprocessing
# .
# Validate 1D Model
# .
# Conservations: Implement 2D equations


log.updateLog()
log.export()
