from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
from classes.Integrator.Integrator import Integrator
from classes.Log.Log import Log
from classes.Postprocessing.Postprocessor import Postprocessor

log = Log("first simulation 1D")

reactor = FixedBedReactor(log, FixedBedReactor.ONE_D, 100)
reactor.setUp()

integrator = Integrator(log, reactor)
integrator.setup(1e-13, 1e-13, 0, 1, 100)
integrator.integrate()

postprocessor = Postprocessor()


#TODO
# add all material properties to FixedBedReactor()
# Integrator                    WORKING -> Implement for 2D
# .
# BUGFIXING EQUATIONS
# .
# Postprocessing                TBD
# .
# Validate 1D Model             TBD
# .
# Create ULM diagramm           TBD
# .
# Conservations: Implement physics as CasADi functions ("D)
#       Pressure drop            TBD
#       Mass conservation        TBD
#       Species Conservation     TBD
#       Energy Conservation      TBD




log.updateLog()
log.export()
