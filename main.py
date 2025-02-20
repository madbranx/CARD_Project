"""
    Project of the Course Computer Aided Reactor Design (CARD) of the
    Institute of Chemical Process Engineering (CVT) at
    Karlsruhe Institute of Technology (KIT)

    Project Timespan 12/2024 - 02/2025

    Authors: Till Kasselmann, Maximilian Brand (Group 1)
    till.kasselmann@student.kit.edu
    maximilian.brand@student.kit.edu
"""
###########################################################################################
import numpy as np
from classes.Postprocessing.Studies import Studies

studies = Studies()

###########################################################################################

""" 1) Discretization  - DONE"""
# studies.discretization_study(
#     "discretization",
#     60,
#     500,
#     1000,
#     [500, 300, 200, 150, 125, 100, 80, 60, 50, 40, 30, 20, 15, 10, 5, 3, 2],
#     3000,
#     250,
#     30,
#     [25, 20, 15, 12, 10, 8, 6, 5, 4, 3, 2],
#     False
# )

# studies.measure_time_discretization(120, 12, 10000, 3000)


""" 2) Validation - DONE 
For validation the radial heat transfer coefficient must be set to provided constant value in FixedBedReactor.py (line 194)"""
# studies.validation(
#     "validation",
#     120,
#     12,
#     3000,
#     300,
#     3000
# )


""" 3) Base Case - DONE """
# studies.base_case_2D(
#     "base_case_2D",
#     120,
#     12,
#     3000,
#     10000
# )


""" 4) Ignition and Extinction Arcs - DONE """
# studies.arcs(
#     "arcs_1D",
#     1,
#     3000,
#     300,
#     np.linspace(300, 550, 25),
#     150,
#     log=False
# )

# studies.arcs(
#     "arcs_2D",
#     2,
#     3000,
#     10000,
#     np.linspace(300, 550, 25),
#     120,
#     12,
#     log=True
# )

# studies.combinedArcs(
#     "arcs_2D_1D",
#     3000,
#     10000,
#     np.linspace(300, 550, 25),
#     120,
#     12,
#     log=False
# )


""" 5) Catalyst Parameter Variations - DONE """
# studies.cat_variation_diameter("cat_variation_diameter",
#                       120,
#                       12,
#                       500,
#                       10000,
#                       np.linspace(0.0014, 0.004, 20), # d_cat
#                       log = False
# )

# studies.cat_variation_pore("cat_variation_pore",
#                       120,
#                       12,
#                       500,
#                       10000,
#                       np.linspace(1e-9, 50e-9, 20),  # d_pores
#                       log = False
# )


""" 6) Ignition Behavior - DONE """
# studies.ignition_behavior(
#     "ignition behaviour",
#     120,
#     12,
#     3000,
#     1000,
#     [25, 150, 200, 250, 325, 450, 750, 3000],
#     0.90,
#     0.70,
#     log=False
# )


###########################################################################################


from classes.FixedBedReactor import FixedBedReactor
import casadi as CasADi

reactor = FixedBedReactor(1, 100, 1)
reactor.setup()

w = np.array([0.5, 0.5, 0., 0.])
T = 500
p = 5e5
u = 1.5
u_center = 1.5
wTpu = [w, T, p, u]
reactorDiameter = 0.02

rad_therm_cond = reactor.calc_effective_radial_thermal_conductivity(wTpu, u_center, 0.0)
print("Radial lambda: ", rad_therm_cond)
alpha_wall = reactor.calc_heatTransferCoefficient_wall(wTpu)
print("Alpha wall: ", alpha_wall)
overall = (1/alpha_wall + reactorDiameter/(2*rad_therm_cond) )**(-1)
print("Overall: ", overall)
