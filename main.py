"""
    Project of the Course Computer Aided Reactor Design (CARD) of the
    Institute of Chemical Process Engineering (CVT) at
    Karlsruhe Institute of Technology (KIT)

    Project Timespan 12/2024 - 02/2025

    Authors: Till Kasselmann, Maximilian Brand (Group 1)
    maximilian.brand@student.kit.edu
    till.kasselmann@student.kit.edu
"""
###########################################################################################

from classes.Studies import Studies
import numpy as np

studies = Studies()

###########################################################################################

#n_axials = [500, 400, 300, 200, 150, 125, 100, 80, 60, 50, 40, 30, 20, 15, 10, 5, 3, 2]
#studies.discretization_study1D(500, 1000, 1000, n_axials)

#n_radials = [25, 20, 15, 12, 10, 8, 6, 4, 3, 2]
#studies.discretization_study_radial(500, 10000, 30, 200, n_radials)

#TODO equi- vs non-equidistant sim time comparison 1D/2D

#studies.validation_1D(300, 3000, 500)

#studies.arcs_1d(150, 4000, 1000, np.linspace(300, 550, 25))

#studies.plot_2D_with_TwallVal(150, 12, 3000, 6000)

#studies.pseudo_2D_vs_1D(20, 5, 3000, 300)

#studies.test_2D_extinction(20, 4, 3000, 8000, 550, 500)

#studies.arcs_2d(150, 12, 3000, 30000, np.linspace(300, 550, 50))

#TODO catalyst variations

#TODO moving Hotspot

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
