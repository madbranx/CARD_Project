from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
from classes.Log.Log import Log

log = Log("first simulation")
reactor = FixedBedReactor(log, FixedBedReactor.ONE_D, 10)
reactor.setUp()


# TODO
# FixedBed: RSQ: define all Parameters, Components (Cat as Component?) & add eps calculation as method
# Conservations: Implement physics as CasADi functions
# FixedBed -> create DAE struct (depending on DIM)


# Integrator
# Postprocessing


log.updateLog()
log.export()
