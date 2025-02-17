from classes.Postprocessing.TOLcmaps import TOLcmaps

import numpy as np
import pandas as pd
import casadi as CasADi
from matplotlib import pyplot as plt
from matplotlib import colors as clr
from scipy.interpolate import interp1d, CubicSpline
from scipy.optimize import curve_fit

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

    def __getValidationData_T_center(self, filename):
        data = pd.read_csv(filename, sep="\t", decimal=",")
        T = [data['z T_center'], data['T_center [K]']]
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



    def plot2D_wTpu_X(self, name, result, timestep):
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

        # Plot Simulation Data Wall & Center
        w_i, T, p, u = result.get_2D_values(timestep)
        z_pos = result.get_z_pos() / self.reactor.reactorLength

        axs.plot(z_pos, T[:, 0], color=colors[1], linestyle="--", label="simulation center")
        axs.plot(z_pos, T[:, -1], color=colors[0], linestyle="--", label="simulation wall")

        # Plot Validation Data Wall
        T_wall_val = self.__getValidationData_T_wall("debugging/Val_data_with_Twall.txt")

        axs.plot(T_wall_val[0], T_wall_val[1], color=colors[0], linestyle="-", label ="validation wall")

        # Plot Validation Data Center

        T_center_val = self.__getValidationData_T_center("debugging/Val_data_with_Twall.txt")

        axs.plot(T_center_val[0], T_center_val[1], color=colors[1], linestyle="-", label ="validation center")

        # Axis
        axs.set_ylabel(r'$T_{\mathregular{wall}} / K$')
        axs.set_xlabel(r'$z/L$')

        axs.legend()

        plt.show()

    def plot_ignitionArc2D(self, results_ignition, results_extinction, T_wall_ign, T_wall_ext, timestep):
        self.setSizes()

        fig, axs = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))


        # Plot Bremer Data
        T_bremer_ignition = np.array([311, 323, 342, 359, 377, 401, 424, 442, 454, 460, 465, 490, 512, 537, 561])
        X_bremer_ignition = np.array([0.00333, 0.00355, 0.00388, 0.00418, 0.00776, 0.00980, 0.0298, 0.0871, 0.201,
                             0.392, 0.930, 0.965, 0.980, 0.985, 0.979])

        T_bremer_extinction = np.array([308, 322, 336, 348, 359, 365, 371, 377, 388, 407, 425, 449, 471, 491, 513, 531, 549])
        X_bremer_extinction = np.array([0.00329, 0.00353, 0.00214, 0.00399, 0.00417, 0.686, 0.826, 0.845, 0.860,
                     0.882, 0.902, 0.931, 0.958, 0.971, 0.985, 0.985, 0.984])

        axs.plot(T_bremer_ignition, X_bremer_ignition, color="black", linestyle="-")
        axs.plot(T_bremer_extinction, X_bremer_extinction, color="black", linestyle="-", linewidth=2)


        # Plot Simulation Data
        X_CO2_ign = []
        for result in results_ignition:
            x_CO2_ign = result.getConversion_2D(timestep, 2)
            X_CO2_ign.append(result.average_trapezoidal(x_CO2_ign[-1, :]))

        X_CO2_ext = []
        for result in results_extinction:
            x_CO2_ext = result.getConversion_2D(timestep, 2)
            X_CO2_ext.append(result.average_trapezoidal(x_CO2_ext[-1, :]))


        axs.plot(T_wall_ign, X_CO2_ign, color=colors[0], linestyle="--", linewidth = 2,  marker='v',
                 markersize=6, markerfacecolor="black",markeredgecolor='black', label = "ignition")
        axs.plot(T_wall_ext, X_CO2_ext, color=colors[1], linestyle="--", linewidth=2, marker='v',
                 markersize=6, markerfacecolor="black", markeredgecolor='black', label = "extinction")

        # Axis
        axs.set_xlabel(r'$T_{\mathregular{wall}} / K$')
        axs.set_ylabel(r'$Conversion CO2 / -$')

        #axs.set_xlim(300, 700)

        axs.legend()
        plt.show()


    def plot1D_vsPseudo1D_vs_val(self, name, result1D, resultP1D, timestep):
        fig, axs = plt.subplots(4, 1, figsize=(4.2, 5.7), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        # Plot Simulation Data 1D
        w_i, T, p, u = result1D.get_values(timestep)
        z_pos = result1D.get_z_pos() / self.reactor.reactorLength

        axs[0].plot(z_pos, w_i[:, 0], color=colors[0], linestyle="--", label="CH4")
        axs[0].plot(z_pos, w_i[:, 1], color=colors[1], linestyle="--", label="H2O")
        axs[0].plot(z_pos, w_i[:, 2], color=colors[2], linestyle="--", label="CO2")
        axs[0].plot(z_pos, w_i[:, 3], color=colors[3], linestyle="--", label="H2")
        axs[1].plot(z_pos, T, color=colors[0], linestyle="--")
        axs[2].plot(z_pos, u, color=colors[0], linestyle="--")
        axs[3].plot(z_pos, p * 1e-5, color=colors[0], linestyle="--", label="1D")

        # Plot Simulation Data P1D
        w_i_2D, T_2D, p_2D, u_2D = resultP1D.get_values(timestep)
        z_pos_2D = resultP1D.get_z_pos() / self.reactor.reactorLength

        axs[0].plot(z_pos_2D, w_i_2D[:, 0], color=colors[0], linestyle=":")
        axs[0].plot(z_pos_2D, w_i_2D[:, 1], color=colors[1], linestyle=":")
        axs[0].plot(z_pos_2D, w_i_2D[:, 2], color=colors[2], linestyle=":")
        axs[0].plot(z_pos_2D, w_i_2D[:, 3], color=colors[3], linestyle=":")
        axs[1].plot(z_pos_2D, T_2D, color=colors[0], linestyle=":")
        axs[2].plot(z_pos_2D, u_2D, color=colors[0], linestyle=":")
        axs[3].plot(z_pos_2D, p_2D * 1e-5, color=colors[0], linestyle=":", label = "P1D")

        # Plot Validation Data
        w_i_val, T_val, u_val, p_val = self.__getValidationData("debugging/Validation_Data.txt")

        axs[0].plot(w_i_val[0][0], w_i_val[0][1], color=colors[0], linestyle="-")
        axs[0].plot(w_i_val[1][0], w_i_val[1][1], color=colors[1], linestyle="-")
        axs[0].plot(w_i_val[2][0], w_i_val[2][1], color=colors[2], linestyle="-")
        axs[0].plot(w_i_val[3][0], w_i_val[3][1], color=colors[3], linestyle="-")
        axs[1].plot(T_val[0], T_val[1], color=colors[0], linestyle="-")
        axs[2].plot(u_val[0], u_val[1], color=colors[0], linestyle="-")
        axs[3].plot(p_val[0], p_val[1], color=colors[0], linestyle="-", label="val")

        axs[3].legend()
        # Axis
        axs[0].set_ylabel(r'$w_{\mathregular{A}}$')
        axs[1].set_ylabel(r'$T \; \mathregular{/K}$')
        axs[2].set_ylabel(r'$u \; \mathregular{/ms^{-1}}$')
        axs[3].set_ylabel(r'$p / bar$')
        axs[3].set_xlabel(r'$z/L$')

        axs[0].legend()

        plt.show()


    def check_sum_wi(self, result, tol):
        for timestep in result.timeDiscretization.get_faces():
            w_i, T, p, u = result.get_2D_values(int(timestep))
            for z in range(len(result.axialDiscretization.get_centroids())):
                for r in range(len(result.radialDiscretization.get_centroids())):
                    summ = 0
                    for i in range(len(w_i)):
                        summ += w_i[i][z][r]
                    if abs(1-summ) >= tol:
                        print("deviation in sum w_i at timestep " + str(timestep) + " at coordinate (" + str(z) + "," +str(r) + "): sum(w_i) = " + str(summ))

        # import pandas as pd
        # # Beispiel-Daten (ersetze dies mit deinen echten Arrays)
        # w_values = np.array(w_i)  # Stoffmengen / Massenanteile
        # T_values = np.array(T)  # Temperatur
        # p_values = np.array(p)  # Druck
        # u_values = np.array(u)  # Geschwindigkeit
        #
        # # Daten als Spalten zusammenfügen
        # data = np.column_stack(T_values)
        #
        # # DataFrame erstellen
        # df = pd.DataFrame(data)
        #
        # # Speichern als CSV
        # df.to_csv("simulation_data.csv", index=False)
        #


    def plot_ignitionArc1D(self, results_ignition, results_extinction, T_walls, time_steps, rev=False):
        self.setSizes()

        fig, axs = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))


        # Plot Bremer Data
        T_bremer_ignition = np.array([311, 323, 342, 359, 377, 401, 424, 442, 454, 460, 465, 490, 512, 537, 561, 596])
        X_bremer_ignition = np.array([0.00333, 0.00355, 0.00388, 0.00418, 0.00776, 0.00980, 0.0298, 0.0871, 0.201,
                             0.392, 0.930, 0.965, 0.980, 0.985, 0.979, 0.970])

        T_bremer_extinction = np.array([308, 322, 336, 348, 359, 365, 371, 377, 388, 407, 425, 449, 471, 491, 513, 531, 549, 573, 595])
        X_bremer_extinction = np.array([0.00329, 0.00353, 0.00214, 0.00399, 0.00417, 0.686, 0.826, 0.845, 0.860,
                     0.882, 0.902, 0.931, 0.958, 0.971, 0.985, 0.985, 0.984, 0.978, 0.971])

        axs.plot(T_bremer_ignition, X_bremer_ignition, color="black", linestyle="-")

        if rev:
            axs.plot(T_bremer_extinction, X_bremer_extinction, color="black", linestyle="-", linewidth=2)

        axs.plot(T_bremer_extinction, X_bremer_extinction, color="black", linestyle="-", linewidth=2)


        # Plot Simulation Data
        X_CO2_ign = []
        for result in results_ignition:
            x_CO2_ign = result.getConversion_1D(time_steps, 2)
            X_CO2_ign.append(x_CO2_ign[-1])

        X_CO2_ext = []
        for result in results_extinction:
            x_CO2_ext = result.getConversion_1D(time_steps, 2)
            X_CO2_ext.append(x_CO2_ext[-1])

        axs.plot(T_walls, X_CO2_ign, color=colors[0], linestyle="--", linewidth = 2,  marker='v', markersize=6, markerfacecolor="black",markeredgecolor='black', label = "ignition")
        axs.plot(np.flip(T_walls), X_CO2_ext, color=colors[1], linestyle="--", linewidth=2, marker='v',
                 markersize=6, markerfacecolor="black", markeredgecolor='black', label = "extinction")

        # Axis
        axs.set_xlabel(r'$T_{\mathregular{wall}} / K$')
        axs.set_ylabel(r'$Conversion CO2 / -$')

        #axs.set_xlim(300, 700)

        axs.legend()
        plt.show()


    def plot_discretizationStudy1D(self, results, result_ref, timesteps):

        ## Calculating RMSEs for every Result
        start, stop = 0, self.reactor.reactorLength
        mean_nRMSEs = []
        max_nRMSEs = []

        # iterating through every timestep in every result
        for result in results:
            centroids_ref = result_ref.reactor.axial_discretization.get_centroids()
            centroids = result.reactor.axial_discretization.get_centroids()

            mean_nRMSEs_timesteps = []
            max_nRMSEs_timesteps = []
            step = int(timesteps / 100)
            for timestep in range(0, timesteps, step):
                mean_nRMSEs_variables = []
                values = result.get_values(timestep)
                ref_values = result_ref.get_values(timestep)
                # Calculating nRMSEs for each Variable (4x w, T, p, u) for one timestep

                ref_values = [ref_values[0][:, 0], ref_values[0][:, 1], ref_values[0][:, 2], ref_values[0][:, 3], ref_values[1][:], ref_values[2][:], ref_values[3][:] ]
                values = [values[0][:, 0], values[0][:, 1], values[0][:, 2], values[0][:, 3], values[1][:], values[2][:], values[3][:] ]

                for i in range(len(ref_values)):
                    values_interp = self.__interpolate_arrays(ref_values[i], centroids_ref, start, stop)
                    ref_values_interp = self.__interpolate_arrays(values[i], centroids, start, stop)
                    mean_nRMSEs_variables.append(self.nRMSE(ref_values_interp,values_interp))


                mean_nRMSEs_timesteps.append(np.mean(mean_nRMSEs_variables))
                max_nRMSEs_timesteps.append(np.max(mean_nRMSEs_variables))
            mean_nRMSEs.append(np.mean(mean_nRMSEs_timesteps))
            max_nRMSEs.append(np.max(max_nRMSEs_timesteps))

        # Plotting
        self.setSizes()

        fig, axs = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        n_axials = []
        for i, result in enumerate(results):
            n_axial = len(result.reactor.axial_discretization.get_centroids())
            n_axials.append(n_axial)
            if i == 0:
                axs.plot(n_axial, mean_nRMSEs[i], marker='v', markersize=10, label='mean', color=colors[0])
                axs.plot(n_axial, max_nRMSEs[i], marker='v', markersize=10, label='max', color=colors[1])

            else:
                axs.plot(n_axial, mean_nRMSEs[i], marker='v', markersize=10, color=colors[0])
                axs.plot(n_axial, max_nRMSEs[i], marker='v', markersize=10, color=colors[1])


        # Create spline fit
        n_axials = np.flip(n_axials)
        spline_mean = CubicSpline(n_axials, np.flip(mean_nRMSEs))
        spline_max = CubicSpline(n_axials, np.flip(max_nRMSEs))

        # Generate smooth x values and compute smooth y values
        x_smooth = np.linspace(min(n_axials), max(n_axials), 100)
        y_smooth_mean = spline_mean(x_smooth)
        y_smooth_max = spline_max(x_smooth)

        plt.plot(x_smooth, y_smooth_mean, linestyle = '--', linewidth = 2, color=colors[0])
        plt.plot(x_smooth, y_smooth_max, linestyle = '--', linewidth = 2, color=colors[1])


        plt.ylabel("normalized RMSE / -")
        plt.xlabel("number of axial volumes/ -")
        plt.legend(loc="upper right")
        plt.title("nRMSE of each variable (4x w_i, T, p, u) \n simulating 500s with tol = 1e-8, reference n_axial = " + str(len(ref_values[0])), fontsize = 15)
        plt.show()

    def nRMSE(self,values_ref, values):
        RMSE = self.RMSE(values_ref, values)
        return RMSE/(np.mean(values))

    def RMSE(self, values_ref, values):
        RMSE = np.sqrt(1/len(values)* np.sum((values- values_ref)** 2))
        return RMSE

    def __interpolate_arrays(self, values, centroids, start, stop):
        values = np.array(values)
        centroids = np.array(centroids)
        resolution = 100
        common_orts = np.linspace(start, stop, resolution)

        interp1 = interp1d(centroids, values, kind='linear', fill_value='extrapolate')

        return np.column_stack((common_orts, interp1(common_orts)))

    def plot_discretizationStudy_radial(self, results, result_ref, timesteps):

        ## Calculating RMSEs for every Result
        start, stop = 0, self.reactor.reactorDiameter/2
        mean_nRMSEs = []
        max_nRMSEs = []


        # iterating through every timestep in every result
        for result in results:
            centroids_ref = result_ref.reactor.radial_discretization.get_centroids()
            centroids = result.reactor.radial_discretization.get_centroids()

            mean_nRMSEs_timesteps = []
            max_nRMSEs_timesteps = []
            step = int(timesteps / 100)
            for timestep in range(0, timesteps, step):
                mean_nRMSEs_axials = []
                max_nRMSEs_axials = []
                for axial in range(len(result_ref.reactor.axial_discretization.get_centroids())):
                    mean_nRMSEs_variables = []
                    values = result.get_2D_values(timestep)
                    ref_values = result_ref.get_2D_values(timestep)
                    # Calculating nRMSEs for each Variable (4x w, T, p, u) for one timestep

                    ref_values = [ref_values[0][0, axial, :], ref_values[0][1, axial, :], ref_values[0][2, axial, :], ref_values[0][3, axial, :], ref_values[1][axial, :], ref_values[2][axial, :], ref_values[3][axial, :] ]
                    values = [values[0][0, axial, :], values[0][1, axial, :], values[0][2, axial, :], values[0][3, axial, :], values[1][axial, :], values[2][axial, :], values[3][axial, :] ]

                    for i in range(len(ref_values)):
                        values_interp = self.__interpolate_arrays(ref_values[i], centroids_ref, start, stop)
                        ref_values_interp = self.__interpolate_arrays(values[i], centroids, start, stop)
                        mean_nRMSEs_variables.append(self.nRMSE(ref_values_interp,values_interp))

                    mean_nRMSEs_axials.append(np.mean(mean_nRMSEs_variables))
                    max_nRMSEs_axials.append(np.max(mean_nRMSEs_variables))
                mean_nRMSEs_timesteps.append(np.mean(mean_nRMSEs_axials))
                max_nRMSEs_timesteps.append(np.max(max_nRMSEs_axials))
            mean_nRMSEs.append(np.mean(mean_nRMSEs_timesteps))
            max_nRMSEs.append(np.max(max_nRMSEs_timesteps))

        # Plotting
        self.setSizes()

        fig, axs = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        n_radials = []
        for i, result in enumerate(results):
            n_radial = len(result.reactor.radial_discretization.get_centroids())
            n_radials.append(n_radial)
            if i == 0:
                axs.plot(n_radial, mean_nRMSEs[i], marker='v', markersize=10, label='mean', color=colors[0])
                axs.plot(n_radial, max_nRMSEs[i], marker='v', markersize=10, label='max', color=colors[1])

            else:
                axs.plot(n_radial, mean_nRMSEs[i], marker='v', markersize=10, color=colors[0])
                axs.plot(n_radial, max_nRMSEs[i], marker='v', markersize=10, color=colors[1])


        # Create spline fit
        n_radials = np.flip(n_radials)
        spline_mean = CubicSpline(n_radials, np.flip(mean_nRMSEs))
        spline_max = CubicSpline(n_radials, np.flip(max_nRMSEs))

        # Generate smooth x values and compute smooth y values
        x_smooth = np.linspace(min(n_radials), max(n_radials), 100)
        y_smooth_mean = spline_mean(x_smooth)
        y_smooth_max = spline_max(x_smooth)

        plt.plot(x_smooth, y_smooth_mean, linestyle = '--', linewidth = 2, color=colors[0])
        plt.plot(x_smooth, y_smooth_max, linestyle = '--', linewidth = 2, color=colors[1])


        plt.ylabel("normalized RMSE / -")
        plt.xlabel("number of radial volumes/ -")
        plt.legend(loc="upper right")
        plt.title("nRMSE of each variable (4x w_i, T, p, u) of each axial row\n simulating 500s with tol = 1e-4, reference n_axial = " + str(len(ref_values[0])), fontsize = 15)
        plt.show()

    def plot2D_T_half(self, name, result, timestep):
        self.setSizes()
        fig, axs = plt.subplots(1, 1, figsize=(20, 5))

        T = result.get_2D_values(timestep)[1]

        z_pos = result.get_z_pos() / self.reactor.reactorLength
        r_pos = 2 * result.get_r_pos() / self.reactor.reactorDiameter

        tol_cmap = TOLcmaps()
        cmp = tol_cmap.get('sunset')

        z_mesh, r_mesh = np.meshgrid(z_pos, r_pos, indexing='ij')

        T_min, T_max = self.get_min_max(T)

        axs.pcolormesh(z_mesh, r_mesh, T, cmap=cmp, vmin=T_min, vmax=T_max,shading='gouraud')

        axs.set_axis_off()  # Hides everything: axes, ticks, and labels
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)  # Removes padding

        plt.savefig("results/T_half_plot.svg", format="svg", bbox_inches="tight")
        plt.savefig("results/T_half_plot.pdf", format="pdf", bbox_inches="tight")
        plt.savefig("results/T_half_plot.eps", format="eps", bbox_inches="tight")
        plt.savefig("results/T_half_plot.png", format="png", dpi=600, bbox_inches="tight")

        plt.show()
