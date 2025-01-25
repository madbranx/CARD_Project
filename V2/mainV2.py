from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator

reactor = FixedBedReactor(1, 100)
integrator = Integrator(reactor)
integrator.setup(1e-13, 1e-13, 0, 10, 100)
integrator.integrate()


#TODO V2
# add collision area for h2o
# Integrator                    WORKING -> Implement for 2D
# .

# .
# Postprocessing                TBD
# .
# Validate 1D Model             TBD
# .
# Create ULM diagramm           TBD
# add Log?
# .
# Conservations: Implement physics as CasADi functions ("D)
#       Pressure drop            TBD
#       Mass conservation        TBD
#       Species Conservation     TBD
#       Energy Conservation      TBD
