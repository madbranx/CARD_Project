from classes.Postprocessing.TOLcmaps import TOLcmaps

import numpy as np
import pandas as pd
import casadi as CasADi
from matplotlib import pyplot as plt
from matplotlib import colors as clr

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

    def __getValidationData_T_wall(self, filename):
        data = pd.read_csv(filename, sep="\t", decimal=",")
        T = [data['z T_wall'], data['T_wall [K]']]
        return T


    def plot_1D_vs_ValidationData(self, name, result, timestep):
        fig, axs = plt.subplots(4, 1, figsize=(4.2, 5.7), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        # Plot Simulation Data
        w_i, T, p, u = result.get_values(timestep)
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

    def plot1D(self, name, result, timestep):

        w_i, T, u, p = result.get_values(timestep)
        z_pos = result.get_z_pos() / self.reactor.reactorLength

        w_f = CasADi.SX.sym('w_f', 4)
        T_f = CasADi.SX.sym('T_f')
        p_f = CasADi.SX.sym('p_f')
        u_f = CasADi.SX.sym('u_f')

        Reynolds = CasADi.Function('f_Re_CasADi', [w_f, T_f, u_f, p_f], [self.reactor.Re(w_f, T_f,u_f, p_f)])

        check = Reynolds(w_i[z_pos, timestep, :], T[z_pos, timestep], p[z_pos, timestep], u[z_pos, timestep])

        fig, ax = plt.subplots(1, 1, figsize=(4.2, 5.7), constrained_layout=True, sharex=True)

        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        ax.plot(z_pos, check, color=colors[0], linestyle="--")

    def setSizes(self):
        # SIZES
        SMALL_SIZE = 18
        MEDIUM_SIZE = 22
        BIGGER_SIZE = 24

        plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



    def plot2D_Temperature(self, name, result, timestep):
        self.setSizes()
        fig, axs = plt.subplots(4, 1, figsize=(7*1.5, 14), constrained_layout=True)

        w_i, T, p, u = result.get_2D_values(timestep)
        x_i = result.getConversion_2D(timestep, 2) # 2 = CO2

        z_pos = result.get_z_pos() / self.reactor.reactorLength
        r_pos = 2*result.get_r_pos() / self.reactor.reactorDiameter
        # Using meshgrid to create coordinate matrices
        z_mesh, r_mesh = np.meshgrid(z_pos, r_pos, indexing='ij')

        # Get color map from TOLcmaps
        tol_cmap = TOLcmaps()
        cmp = tol_cmap.get('sunset')

        # Dynamically determine min/max for colormaps with rounding
        u_min, u_max = self.get_min_max(u)
        p_min, p_max = self.get_min_max(p * 1e-5)  # Convert to bar before min/max
        T_min, T_max = self.get_min_max(T)
        w_i0_min, w_i0_max = self.get_min_max(w_i[0])
        x_i_min, x_i_max = self.get_min_max(x_i)

        def plot_variable(ax, data, cmap, min_val, max_val, label, fmt_color='black'):
            ax.pcolormesh(z_mesh, r_mesh, data, cmap=cmap, vmin=min_val, vmax=max_val)
            plot = ax.pcolormesh(z_mesh, -r_mesh, data, cmap=cmap, vmin=min_val, vmax=max_val)
            levels = self.get_contour_levels(min_val, max_val, data)
            c_lines = ax.contour(z_mesh, r_mesh, data, levels=levels, colors=fmt_color, linewidths=0.8)
            c_lines = ax.contour(z_mesh, -r_mesh, data, levels=levels, colors=fmt_color, linewidths=0.8)
            ax.clabel(c_lines, fmt=self.get_contour_label_format(levels), fontsize=8)

            # Draw grid lines at dataset points
            for z in z_pos:
                ax.axvline(z, color="white", linestyle="--", linewidth=0.5, alpha=0.5)
            for r in r_pos:
                ax.axhline(r, color="white", linestyle="--", linewidth=0.5, alpha=0.5)
                ax.axhline(-r, color="white", linestyle="--", linewidth=0.5, alpha=0.5)

            fig.colorbar(plot, label=label, ax=ax)


        # Plot velocity
        plot_variable(axs[0], u, cmp, u_min, u_max, 'Velocity / m/s')

        # Plot pressure
        plot_variable(axs[1], p * 1e-5, cmp, p_min, p_max, 'Pressure / bar')

        # Plot temperature (use white contour lines for better visibility)
        plot_variable(axs[2], T, 'inferno', T_min, T_max, 'Temperature / K', fmt_color='white')

        # Plot Conversion
        plot_variable(axs[3], x_i, 'Oranges', x_i_min, x_i_max, 'Conversion CO2 / -')

        plt.show()

        X_i = result.average_trapezoidal(x_i[-1, :])
        print(f"X CO2 (integrated) = {X_i}")

    def get_min_max(self, data):
        raw_min, raw_max = np.nanmin(data), np.nanmax(data)  # Ignore NaNs
        rounded_min = self.round_sensefully(raw_min, 'down')
        rounded_max = self.round_sensefully(raw_max, 'up')
        return rounded_min, rounded_max

    def round_sensefully(self,value, round_type='down'):
        if value == 0:
            return 0  # Avoid log(0) errors

        magnitude = 10 ** np.floor(np.log10(abs(value)))  # Get order of magnitude
        if round_type == 'down':
            return np.floor(value / magnitude) * magnitude
        else:
            return np.ceil(value / magnitude) * magnitude

    def get_contour_levels(self, min_val, max_val, data):
        percentiles = np.percentile(data, [10, 25, 50, 75, 90])
        levels = np.unique(np.concatenate(([min_val], percentiles, [max_val])))
        return levels

    def get_contour_label_format(self, levels):
        max_val = max(abs(levels.min()), abs(levels.max()))

        if max_val >= 100:  # Large values
            return "%.0f"
        elif max_val >= 10:
            return "%.1f"
        elif max_val >= 1:  # Normal values
            return "%.2f"
        else:  # Small values
            return "%.3f"

    def plot_Twall_vs_validation(self, name, result, timestep):
        fig, axs = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        # Plot Simulation Data
        w_i, T, p, u = result.get_2D_values(timestep)
        z_pos = result.get_z_pos() / self.reactor.reactorLength

        axs.plot(z_pos, T[:, 0], color=colors[0], linestyle="--", label="simulation")

        # Plot Validation Data
        T_wall_val = self.__getValidationData_T_wall("debugging/Val_data_with_Twall.txt")

        axs.plot(T_wall_val[0], T_wall_val[1], color=colors[0], linestyle="-", label ="validation values")


        # Axis
        axs.set_ylabel(r'$T_{\mathregular{wall}} / K$')
        axs.set_xlabel(r'$z/L$')

        axs.legend()

        plt.show()

    def plot_ignitionArc(self, results_ignition, results_extinction, T_walls, timestep):
        self.setSizes()

        X_CO2_ign = []
        for result in results_ignition:
            x_CO2_ign = result.getConversion_2D(timestep, 2)
            X_CO2_ign.append(result.average_trapezoidal(x_CO2_ign[-1, :]))

        X_CO2_ext = []
        for result in results_extinction:
            x_CO2_ext = result.getConversion_2D(timestep, 2)
            X_CO2_ext.append(result.average_trapezoidal(x_CO2_ext[-1, :]))

        fig, axs = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        axs.plot(T_walls, X_CO2_ign, color=colors[0], linestyle="-", linewidth = 2,  marker='v', markersize=6, markerfacecolor="black",markeredgecolor='black', label = "ignition")
        axs.plot(T_walls, X_CO2_ext, color=colors[1], linestyle="-", linewidth=2, marker='v',
                 markersize=6, markerfacecolor="black", markeredgecolor='black', label = "extinction")

        # Axis
        axs.set_xlabel(r'$T_{\mathregular{wall}} / K$')
        axs.set_ylabel(r'$Conversion CO2 / -$')

        axs.set_xlim(300, 700)

        axs.legend()
        plt.show()



