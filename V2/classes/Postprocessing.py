














# # TODO below here only for first testing!
#
#
# w_i_in = self.reactor.w_i_in
# T_in = self.reactor.T_in
# u_in = self.reactor.u_in
# p_in = self.reactor.p_in
#
# MassFluxDev = np.empty(shape=(n_spatial, t_steps))
# mdot_0 = (u_in * self.reactor.rho_fl(w_i_in, T_in, p_in)).__float__()
# for t in range(t_steps):
#     for z in range(n_spatial):
#         MassFluxDev[z, t] = abs(mdot_0 - (u_res[z, t] * self.reactor.rho_fl(w_i_res[z, t, :].T, T_res[z, t],
#                                                                             p_res[z, t]).__float__())) / mdot_0 * 100
#         # print(MassFluxDev[z, t])
# print("Maximal mass flux deviation: ", np.max(MassFluxDev))
#
# w_f = CasADi.SX.sym('w_f', 4)
# T_f = CasADi.SX.sym('T_f')
# p_f = CasADi.SX.sym('p_f')
# u_f = CasADi.SX.sym('u_f')
#
# # ETA
# f_eta_casADI = CasADi.Function('f_eta_casADi', [w_f, T_f, p_f], [self.reactor.effFactor(w_f, T_f, p_f)])
# eta = np.empty(shape=(n_spatial, t_steps))
#
# # check function
# f_check_casADI = CasADi.Function('f_check_casADI', [w_f, T_f, p_f, u_f], [])
# for t in range(t_steps):
#     for z in range(n_spatial):
#         if z == 0:
#             pass
#         eta[z, t] = f_eta_casADI(w_i_res[z, t, :], T_res[z, t], p_res[z, t])
#         # print(eta[z, t])
#         check = f_check_casADI(w_i_res[z, t, :], T_res[z, t], p_res[z, t], u_res[z, t])
#         # print("w_i", w_i_res[z, t, :])
#         # print("conv heat flux", self.reactor.axialMassFlow(T_res[z, t], w_i_res[z, t, :], w_i_res[z-1, t, :], u_res[z, t], p_res[z, t], 2))
#         # print("rate eq.", check)
#
# # Plot results
# # Generate plot
# fig, axs = plt.subplots(5, 1, figsize=(4.2, 6.7), constrained_layout=True, sharex=True)
# # Define colors
# colors = plt.cm.Dark2(np.linspace(0, 1, 8))
# # Plot
# t_step = 1000
#
# axs[0].plot(w_i_res[:, t_step, 0], color=colors[0])
# axs[0].plot(w_i_res[:, t_step, 1], color=colors[1])
# axs[0].plot(w_i_res[:, t_step, 2], color=colors[2])
# axs[0].plot(w_i_res[:, t_step, 3], color=colors[3])
# axs[1].plot(T_res[:, t_step], color=colors[0])
# axs[1].set_ylim([300, 1100])
# axs[2].plot(u_res[:, t_step], color=colors[0])
# axs[3].plot(p_res[:, t_step] * 1e-5, color=colors[0])
#
# axs[4].plot(eta[:, t_step], color=colors[0])
# axs[4].set_ylim([0, 1.1])
# # Axis
# axs[0].set_ylabel(r'$w_{\mathregular{A}}$')
# axs[1].set_ylabel(r'$T \; \mathregular{/K}$')
# axs[2].set_ylabel(r'$u \; \mathregular{/ms^{-1}}$')
# axs[3].set_ylabel(r'$p / bar$')
# axs[4].set_ylabel(r'$eff Factor / -$')
# axs[4].set_xlabel(r'$z/L$')
# plt.show()