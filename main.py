from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
from classes.Log.Log import Log

log = Log("first simulation")
reactor = FixedBedReactor(log, FixedBedReactor.ONE_D, 1)
reactor.setUp()


#TODO
# add all parameters and material properties to FixedBedReactor()
# Conservations: Implement physics as CasADi functions (1D)
#       Pressure drop            DONE -> make ergun class ... (same structure as Species Conservation)
#       Mass conservation        DONE -> apply structure here too
#       Species Conservation     DONE
#       Energy Conservation      TBD
# .
# .
# Integrator
# Postprocessing
# .
# Validate 1D Model
# .
# implement 2D equations


log.updateLog()
log.export()
