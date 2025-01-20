from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
from classes.Log.Log import Log

log = Log("first simulation")
reactor = FixedBedReactor(log, FixedBedReactor.ONE_D, 10)


# TODO
# Log: add log messages/warnings/errors to classes                                                                  DONE
# FixedBed: RSQ: define all Parameters, Components (Cat as Component?) & add eps calculation as method
# FixedBed -> Conservations: add disc, dim, RSC as attributes                                                       DONE
# Conservations: Implement physics as CasADi functions
# Discretization: add array functionality                                                                           DONE
# FixedBed -> create DAE struct (depending on DIM)

# Maybe integrate DIM, DISC, LOG into RSC?

# Integrator
# Postprocessing


log.updateLog()
log.export()
