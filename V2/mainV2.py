from classes.FixedBedReactor import FixedBedReactor
from classes.Integrator import Integrator

reactor = FixedBedReactor(1, 5)
integrator = Integrator(reactor)
integrator.setup(1e-13, 1e-13, 0, 1, 100)
integrator.integrate()


#TODO V2
# add collision area for h2o
# ADD INPUT / OUTPUT UNITS FOR ALL FUNCTIONS as comments
# .
# Postprocessing                TBD
# .
# Validate 1D Model             TBD
# .
# Create ULM diagram           TBD
# add Log?
# .
# Conservations: Implement physics as CasADi functions ("D)
#       Pressure drop            TBD
#       Mass conservation        TBD
#       Species Conservation     TBD
#       Energy Conservation      TBD
