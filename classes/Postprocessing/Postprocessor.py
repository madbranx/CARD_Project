class Postprocessor:

    def __init__(self):
        pass



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