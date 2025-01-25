from V1.classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
from V1.classes.Integrator.Integrator import Integrator
from V1.classes.Log.Log import Log
from V1.classes.Postprocessing.Postprocessor import Postprocessor

log = Log("first simulation 1D")

reactor = FixedBedReactor(log, FixedBedReactor.ONE_D, 20)
reactor.setUp()

integrator = Integrator(log, reactor)
integrator.setup(1e-13, 1e-13, 0, 1, 100)
integrator.integrate()

postprocessor = Postprocessor()


#TODO
# add collision area for h2o
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
