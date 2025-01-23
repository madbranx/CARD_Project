from classes.FixedBedReactor.FixedBedReactor import FixedBedReactor
from classes.Integrator.Integrator import Integrator
from classes.Log.Log import Log
from classes.Postprocessing.Postprocessor import Postprocessor

log = Log("first simulation")

reactor = FixedBedReactor(log, FixedBedReactor.ONE_D, 2)
reactor.setUp()

integrator = Integrator(reactor)
integrator.integrate()

postprocessor = Postprocessor()


#TODO
# add all material properties to FixedBedReactor()
# Integrator
# Postprocessing
# .
# Validate 1D Model
# .
# Conservations: Implement physics as CasADi functions ("D)
#       Pressure drop            TBD
#       Mass conservation        TBD
#       Species Conservation     TBD
#       Energy Conservation      TBD

# # Define integrator method
# options = {'abstol': 1e-13, 'reltol': 1e-13}
# I = casADi.integrator('I', 'idas', dae, t_0, t_i, options)
# # Integrate
# x_0 = np.zeros(n_axial * (len(nu) + 1))
# x_0[0:n_axial] = w_i_0[0]
# x_0[n_axial:2 * n_axial] = w_i_0[1]
# x_0[2 * n_axial:] = T_0
# z_0= np.zeros(n_axial * 2)
# z_0[0:n_axial] = u_in
# z_0[n_axial:] = p_in
# res_dic = I(x0=x_0, z0=z_0)
#
# # Reshape results
# res_x = res_dic['xf'].full()
# w_i_res = np.empty(shape=(n_axial, t_steps, len(nu)))
# T_res = np.empty(shape=(n_axial, t_steps))
# for t in range(t_steps):
#     T_res[:, t] = res_x[len(nu) * n_axial:, t]
#     for comp in range(len(nu)):
#         w_i_res[:, t, comp] = res_x[comp * n_axial: n_axial * (comp + 1), t]
# ae_res = res_dic['zf'].full()
# u_res = ae_res[:n_axial, :]
# p_res = ae_res[n_axial:, :]



log.updateLog()
log.export()
