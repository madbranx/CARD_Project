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
from classes.Studies import Studies

studies = Studies()


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
#     True
# )

#TODO equi- vs non-equidistant sim time comparison 1D/2D


""" 2) Validation - DONE """
# studies.validation(
#     "validation",
#     150,
#     12,
#     3000,
#     300,
#     15000
# )


""" 3) Base Case - DONE """
# studies.base_case_2D(
#     "base_case_2D",
#     200,
#     20,
#     1500,
#     12000
# )


""" 4) Ignition and Extinction Arcs - DONE """
# studies.arcs(
#     "arcs_1D",
#     1,
#     450,
#     100,
#     np.linspace(300, 550, 15),
#     100,
#     log=False
# )

# studies.arcs(
#     "arcs_2D",
#     2,
#     3000,
#     10000,
#     np.linspace(300, 550, 25),
#     80,
#     6,
#     log=True
# )

# studies.combinedArcs(
#     "arcs_2D_1D",
#     100,
#     100,
#     np.linspace(450, 550, 1),
#     20,
#     4,
#     log=False)


""" 5) Catalyst Parameter Variations - DONE """
# studies.cat_variation_diameter("cat_variation_diameter",
#                       20,
#                       4,
#                       500,
#                       2000,
#                       np.linspace(0.0014, 0.004, 20), # d_cat
#                       log = False
# )

# studies.cat_variation_pore("cat_variation_pore",
#                       30,
#                       5,
#                       500,
#                       2000,
#                       np.linspace(1e-9, 50e-9, 20),  # d_pores
#                       log = False
# )


""" 6) Ignition Behavior """
#TODO moving Hotspot





###########################################################################################

#studies.plot_2D_with_TwallVal(150, 12, 3000, 6000)

#studies.pseudo_2D_vs_1D(20, 5, 3000, 300)

#studies.test_2D_extinction(20, 4, 3000, 8000, 550, 500)

###########################################################################################

## TESTING FUNCTIONS

# import casadi as CasADi
# T = 50s0
# p = 5e5
# u = 1ss
#
# w1 = CasADi.SX([0.25, 0.25, 0.25, 0.25])
# w2 = CasADi.SX([0, 0, 0.800001, 0.199999])
#
# wTpu1 = [w1, T, p, u]
# wTpu2 = [w2, T, p, u]
#
# print("j dispersion = ", reactor.calc_j_dispersion(0.1, 0.1, 0.1, wTpu1, wTpu2, 2))
