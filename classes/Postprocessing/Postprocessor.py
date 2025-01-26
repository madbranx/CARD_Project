import numpy as np
import pandas as pd
import casadi as CasADi
from matplotlib import pyplot as plt

class Postprocessor:
    def __init__(self, reactor, exportLocation):
        self.reactor = reactor
        self.exportLocation = exportLocation

    def __getValidationData(self, filename):
        data = pd.read_csv(filename, sep="\t", decimal=",")

        w_CH4 = [data['z x_CH4'], data['w_CH4']]
        w_H2O = [data['z x_H2O'], data['w_H2O']]
        w_CO2 = [data['z x_CO2'], data['w_CO2']]
        w_H2 = [data['z x_H2'], data['w_H2']]

        w_i = [w_CH4, w_H2O, w_CO2, w_H2]
        T = [data['z T_center'], data['T_center [K]']]

        p = [data['z p'], data['p [bar]']]
        u = [data['z u_center'], data['u_center [ms]']]

        return w_i, T, u, p


    def plot_1D_vs_ValidationData(self, name, result, timestep):
        fig, axs = plt.subplots(4, 1, figsize=(4.2, 5.7), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        # Plot Simulation Data
        w_i, T, u, p = result.get_values(timestep)
        z_pos = result.get_z_pos()/self.reactor.reactorLength

        axs[0].plot(z_pos, w_i[:, 0], color=colors[0], linestyle="--", label="CH4")
        axs[0].plot(z_pos, w_i[:, 1], color=colors[1], linestyle="--", label="H2O")
        axs[0].plot(z_pos, w_i[:, 2], color=colors[2], linestyle="--", label="CO2")
        axs[0].plot(z_pos, w_i[:, 3], color=colors[3], linestyle="--", label="H2")
        axs[1].plot(z_pos, T, color=colors[0], linestyle="--")
        axs[2].plot(z_pos, u, color=colors[0], linestyle="--")
        axs[3].plot(z_pos, p * 1e-5, color=colors[0], linestyle="--")

        # Plot Validation Data
        w_i_val, T_val, u_val, p_val = self.__getValidationData("debugging/Validation_Data.txt")

        axs[0].plot(w_i_val[0][0], w_i_val[0][1], color=colors[0], linestyle="-")
        axs[0].plot(w_i_val[1][0], w_i_val[1][1], color=colors[1], linestyle="-")
        axs[0].plot(w_i_val[2][0], w_i_val[2][1], color=colors[2], linestyle="-")
        axs[0].plot(w_i_val[3][0], w_i_val[3][1], color=colors[3], linestyle="-")
        axs[1].plot(T_val[0], T_val[1], color=colors[0], linestyle="-")
        axs[2].plot(u_val[0], u_val[1], color=colors[0], linestyle="-")
        axs[3].plot(p_val[0], p_val[1], color=colors[0], linestyle="-")

        # Axis
        axs[0].set_ylabel(r'$w_{\mathregular{A}}$')
        axs[1].set_ylabel(r'$T \; \mathregular{/K}$')
        axs[2].set_ylabel(r'$u \; \mathregular{/ms^{-1}}$')
        axs[3].set_ylabel(r'$p / bar$')
        axs[3].set_xlabel(r'$z/L$')

        axs[0].legend()

        plt.show()

    def plot1D(self, name, result):
        pass



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
        #         check = f_check_casADI(w_i_res[z, t, :], T_res[z, t], p_res[z, t], u_res[z, t])
        #
        #







