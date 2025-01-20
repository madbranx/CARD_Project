from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
from classes.Log.Log import Log

log = Log("log")
reactor = FixedBedReactor(log, FixedBedReactor.TWO_D, 50, 5)


# TODO
# Log: add log messages/warnings/errors to classes
# FixedBed: RSQ: define all Parameters, Components (Cat as Component?)
# FixedBed -> Conservations add disc, dim, RSC as attributes
# Conservations: Implement physics as CasADi functions
# Discretization: add array/function functionality
# FixedBed -> create DAE struct (depending on DIM)

# Integrator
# Postprocessing
