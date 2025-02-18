import numpy as np
import pandas as pd
import casadi as CasADi
import math
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import interp1d, griddata
from pathlib import Path
from classes.Postprocessing.TOLcmaps import TOLcmaps


class Postprocessor:
    def __init__(self, exportLocation, use_latex = False):
        self.exportLocation = exportLocation
        self.useLatex = use_latex #if True needs MikTex installed on PC!!!

    ''' '################################ FORMAT AND DATA EXTRACTION METHODS ################################'''

    def __getValidationData(self, filename, variable):
        data = pd.read_csv(filename, sep="\t", decimal=",")

        if variable == "T":
            return np.array([data['z T_center'], data['T_center [K]']])
        elif variable == "p":
            return np.array([data['z p'], data['p [bar]']])
        elif variable == "u":
            return np.array([data['z u_center'], data['u_center [ms]']])
        elif variable == "w":
            w_CH4 = [data['z x_CH4'], data['w_CH4']]
            w_H2O = [data['z x_H2O'], data['w_H2O']]
            w_CO2 = [data['z x_CO2'], data['w_CO2']]
            w_H2 = [data['z x_H2'], data['w_H2']]
            return np.array([w_CH4, w_H2O, w_CO2, w_H2])
        elif variable == "T_wall":
            return  np.array([data['z T_wall'], data['T_wall [K]']])

    def __format_plot(self):
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

        # set font type
        if self.useLatex is True:
            plt.rcParams.update({
                                    "pgf.texsystem": "pdflatex",
                                    "font.family": "serif",
                                    "text.usetex": True,
                                    "pgf.rcfonts": False,
                                })
            plt.rcParams["font.serif"] = ["Computer Modern Roman", "Times", "Georgia"]
        else:
            plt.rcParams["font.family"] = "serif"
            plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif", "Georgia"]

    def __get_colors(self):
        #blue purple teal yellow
        return ['#004488', '#AA3377', '#009988', '#DDAA33']

    def __export(self, foldername, plotname, plot):
        # Create directories if not exists
        Path(self.exportLocation).mkdir(parents=True, exist_ok=True)
        folder_path = Path(self.exportLocation) / foldername
        folder_path.mkdir(parents=True, exist_ok=True)

        exportname = folder_path / plotname

        try:
            plot.savefig(f"{exportname}.tiff", format="tiff", bbox_inches="tight")
            plot.savefig(f"{exportname}.svg", format="svg", bbox_inches="tight")
            plot.savefig(f"{exportname}.pdf", format="pdf", bbox_inches="tight", transparent=True)
            plot.savefig(f"{exportname}.eps", format="eps", bbox_inches="tight", transparent=True)
            plot.savefig(f"{exportname}.png", format="png", dpi=600, bbox_inches="tight")
        except Exception as e:
            print(f"Error saving plot: {e}")

    '''##################################### PLOT METHODS  - Validation #####################################'''

    def plot_1D_and_2D_vs_ValidationData(self, foldername, result_1D, result_2D, timestep_1D, timestep_2D):
        colors = self.__get_colors()
        self.__format_plot()
        linewidth = 2.5

        # get values of center cell [:, 0] when results are 2D
        w_i, T, p, u = result_2D.get_2D_values(timestep_2D)
        w_i_2D = [w_i[comp, :, 0] for comp in range(len(result_2D.reactor.components))]
        T_2D, p_2D, u_2D = T[:, 0], p[:, 0]*1e-5, u[:, 0]

        w_i, T, p, u = result_1D.get_values(timestep_1D)

        p = p * 1e-5
        z_pos_1D = result_1D.get_z_pos() / result_1D.reactor.reactorLength
        z_pos_2D = result_2D.get_z_pos() / result_2D.reactor.reactorLength

        # Create side-by-side subplots
        fig, (axs1, axs2) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True, sharex=True)

        ## 1) w_i plot
        for i in range(len(result_2D.reactor.components)):
            axs1.plot(z_pos_1D, w_i[:, i], color=colors[i], linestyle=":", linewidth = linewidth)
            axs1.plot(z_pos_2D, w_i_2D[i], color=colors[i], linestyle="--", linewidth = linewidth)

        # Plot Validation Data
        w_i_val = self.__getValidationData("debugging/Validation_Data.txt", "w")
        for i, label in enumerate(["CH4", "H2O", "CO2", "H2"]):
            axs1.plot(w_i_val[i][0], w_i_val[i][1], color=colors[i], linestyle="-", linewidth = linewidth, label=label)

        axs1.set_ylabel(r'$w_\mathrm{{i}} / -$')
        axs1.set_xlabel(r'$z/L_\mathrm{{reactor}}$')
        axs1.set_yticks([0, 0.5, 1])
        axs1.set_xticks([0, 0.5, 1])
        axs1.set_ylim([0, 1])
        axs1.set_xlim([0, 1])
        axs1.legend(loc="upper right", ncol=2)

        ## 2) T, p, u plot
        T_val = self.__getValidationData("debugging/Validation_Data.txt", "T")
        u_val = self.__getValidationData("debugging/Validation_Data.txt", "u")
        p_val = self.__getValidationData("debugging/Validation_Data.txt", "p")

        # Get normalization values 1D
        T_max = np.max([np.max(T_val[1]), np.max(T)])
        u_max = np.max([np.max(u_val[1]), np.max(u)])
        p_max = np.max([np.max(p_val[1]), np.max(p)])

        # Plot Simulation Data 1D
        axs2.plot(z_pos_1D, T / T_max, color=colors[0], linestyle=":", linewidth = linewidth)
        axs2.plot(z_pos_1D, u / u_max, color=colors[1], linestyle=":", linewidth = linewidth)
        axs2.plot(z_pos_1D, p / p_max, color=colors[2], linestyle=":", linewidth = linewidth)

        # Get normalization values 2D
        T_max_2D = np.max([np.max(T_val[1]), np.max(T_2D)])
        u_max_2D = np.max([np.max(u_val[1]), np.max(u_2D)])
        p_max_2D = np.max([np.max(p_val[1]), np.max(p_2D)])

        # Plot Simulation Data 2D
        axs2.plot(z_pos_2D, T / T_max_2D, color=colors[0], linestyle="--", linewidth = linewidth)
        axs2.plot(z_pos_2D, u / u_max_2D, color=colors[1], linestyle="--", linewidth = linewidth)
        axs2.plot(z_pos_2D, p / p_max_2D, color=colors[2], linestyle="--", linewidth = linewidth)


        # Plot Validation Data
        axs2.plot(T_val[0], T_val[1] / T_max, color=colors[0], linestyle="-", linewidth = linewidth, label=r"$Temperature$")
        axs2.plot(u_val[0], u_val[1] / u_max, color=colors[1], linestyle="-", linewidth = linewidth, label=r"$Velocity$")
        axs2.plot(p_val[0], p_val[1] / p_max, color=colors[2], linestyle="-", linewidth = linewidth, label=r"$Pressure$")

        axs2.set_ylabel(r'normalized value / -')
        axs2.set_xlabel(r'$z/L_{\mathrm{reactor}}$')
        axs2.set_yticks([0, 0.5, 1])
        axs2.set_xticks([0, 0.5, 1])
        axs2.set_ylim([0, 1])
        axs2.set_xlim([0, 1])
        axs2.legend(loc="lower right")

        # Third legend (for explanation of line styles)
        line_styles = [plt.Line2D([0], [0], color="black", linestyle="dotted", lw=linewidth),
                       plt.Line2D([0], [0], color="black", linestyle="dashed", lw=linewidth),
                       plt.Line2D([0], [0], color="black", linestyle="solid", lw=linewidth)]
        labels = ["1D Simulation", "2D Simulation", "Validation values"]

        fig.legend(handles=line_styles, labels=labels, loc="upper center",
                   bbox_to_anchor=(0.5, 1.1), ncol=3)

        #plt.tight_layout()

        # Add Title
        # if self.useLatex:
        #     fig.suptitle(r"\textbf{Comparison of the Simulation with the Validation Data\n}", y=1.125)
        # else:
        #     fig.suptitle("Comparison of the Simulation with Validation Data\n", y=1.125)

        # Export the combined figure
        self.__export(foldername, "sim_validation_normalized", plt)
        plt.show()

    '''##################################### PLOT METHODS  - Base Case #####################################'''

    def plot2D_wTpu_X(self, foldername, result, timestep, conversion_comp):
        self.__format_plot()

        fig, axs = plt.subplots(4, 1, figsize=(8 * 1.6, 8), constrained_layout=True)

        w_i, T, p, u = result.get_2D_values(timestep)
        p = p * 1e-5
        x_i = result.getConversion_2D(timestep, conversion_comp)

        z_pos = result.get_z_pos() / result.reactor.reactorLength
        r_pos = 2 * result.get_r_pos() / result.reactor.reactorDiameter
        z_mesh, r_mesh = np.meshgrid(z_pos, r_pos, indexing='ij')

        # Nested Function for Plotting
        def plot_variable(ax, data, cmap, min_val, max_val, title_text, contour_color='black', rev_contour=False, ticks = None):
            plot = ax.pcolormesh(z_mesh, r_mesh, data, cmap=cmap, vmin=min_val, vmax=max_val)#, shading='gouraud')
            plot = ax.pcolormesh(z_mesh, -r_mesh, data, cmap=cmap, vmin=min_val, vmax=max_val)#, shading='gouraud')

            # colorbar
            cbar = fig.colorbar(plot, ax=ax, aspect=7, shrink=1, pad=0.01, location="right")
            if ticks is not None:
                cbar.set_ticks(ticks)

            # Add titles as texts
            ax.text(0.5, 1.05, title_text, ha='center', va='bottom', fontsize=20, transform=ax.transAxes)

            # Combine r_mesh and -r_mesh and data
            combined_r_pos = np.concatenate([r_mesh, -r_mesh], axis=0)
            combined_data = np.concatenate([data, data], axis=0)
            combined_z_pos = np.concatenate([z_mesh, z_mesh], axis=0)

            # draw contour lines

            # Define the fine grid
            z_pos_fine = np.linspace(combined_z_pos.min(), combined_z_pos.max(), 1000)
            r_pos_fine = np.linspace(combined_r_pos.min(), combined_r_pos.max(), 500)
            z_mesh_fine, r_mesh_fine = np.meshgrid(z_pos_fine, r_pos_fine)

            # Interpolate the data on the finer grid
            points = np.column_stack((combined_z_pos.flatten(), combined_r_pos.flatten()))
            values = combined_data.flatten()
            data_interp = griddata(points, values, (z_mesh_fine, r_mesh_fine), method='cubic')

            levels = self.__get_contour_levels(np.min(data_interp), np.max(data_interp), data_interp, rev_contour)
            c_lines = ax.contour(z_mesh_fine, r_mesh_fine, data_interp, levels=levels, colors=contour_color, linewidths=0.8, antialiased=True)
            ax.clabel(c_lines, fmt=self.__get_contour_label_format(levels), fontsize=12)


        # Get color map from TOLcmaps
        tol_cmap = TOLcmaps()

        # Dynamically determine min/max for colormaps with rounding
        u_min, u_max = self.__get_min_max(u)
        p_min, p_max = self.__get_min_max(p)
        T_min, T_max = self.__get_min_max(T)
        x_i_min, x_i_max = self.__get_min_max(x_i)

        # Plot velocity
        plot_variable(axs[0], u, tol_cmap.get('BuRd'), u_min, u_max, r"u / m/s")

        # Plot pressure
        plot_variable(axs[1], p, tol_cmap.get('BuRd'), p_min, p_max, r'p / bar')

        # Plot temperature (use white contour lines for better visibility)
        plot_variable(axs[2], T, tol_cmap.get('sunset'), T_min, T_max, r'T / K', ticks = [300, 600, 900])

        # Plot Conversion
        comps = [r"CH_4", r"H_20", r"CO_2", r"H_2"]
        plot_variable(axs[3], x_i, "Greys_r", x_i_min, x_i_max, r"$X_{ " + comps[conversion_comp] + "}$ / -", '#AA3377', True, ticks=[0, 0.5, 1])

        # Modify Y-axis
        for ax in axs:
            ax.set_yticks([-1, 0, 1])
            ax.set_yticklabels(['R',"0", 'R'])  # Replace -1 and 1 with "R"
            ax.set_xticks([])
        axs[3].set_xticks([0, 0.5, 1])  # Set x-ticks to 0, 0.5, 1
        axs[3].set_xticklabels([0, 0.5, 1])
        axs[3].tick_params(axis='x', pad=10)

        # Set x-axis label on the last subplot
        axs[-1].set_xlabel("relative reactor length")

        fig.text(-0.035, 0.5, "radial position", va='center', rotation=90, fontsize=22)

        self.__export(foldername, "XTpu_plot", plt)
        plt.show()

        x_i = result.getConversion_2D(timestep, conversion_comp)
        X_i = result.average_trapezoidal(x_i[-1, :])
        print("X " + comps[conversion_comp] + f" (integrated) at outlet = {X_i}")

    ####### additional methods for base case plot ########

    def __get_min_max(self, data):
        raw_min, raw_max = np.nanmin(data), np.nanmax(data)  # Ignore NaNs
        rounded_min = self.__round_sensefully(raw_min, 'down')
        rounded_max = self.__round_sensefully(raw_max, 'up')
        return rounded_min, rounded_max

    def __round_sensefully(self,value, round_type='down'):
        if value == 0:
            return 0  # Avoid log(0) errors

        magnitude = 10 ** np.floor(np.log10(abs(value)))  # Get order of magnitude
        if round_type == 'down':
            return np.floor(value / magnitude) * magnitude
        else:
            return np.ceil(value / magnitude) * magnitude

    def __get_contour_levels(self, min_val, max_val, data, rev=False):
        if rev is True:
            percentiles = np.percentile(data, [15, 40, 70])
        else:
            percentiles = np.percentile(data, [5, 75, 95])
        levels = np.unique(np.concatenate(([min_val], percentiles, [max_val])))
        return levels

    def __get_contour_label_format(self, levels):
        max_val = max(abs(levels.min()), abs(levels.max()))

        if max_val >= 100:  # Large values
            return "%.0f"
        elif max_val >= 10:
            return "%.1f"
        elif max_val >= 1:  # Normal values
            return "%.2f"
        else:  # Small values
            return "%.3f"

    '''################################ PLOT METHODS  - Discretization Study ################################'''

    def plot_disdiscretizationStudy(self, foldername, result_axial_ref, results_axial_ed, results_axial_ned, result_radial_ref, results_radial_ed, results_radial_ned, timesteps_axial, timesteps_radial):
        # calculating nRMSE values
        ax_ed_mean, ax_ed_max = self.__getRMSE_1D(results_axial_ed, result_axial_ref, timesteps_axial)
        ax_ned_mean, ax_ned_max = self.__getRMSE_1D(results_axial_ned, result_axial_ref, timesteps_axial)
        rad_ed_mean, rad_ed_max = self.__getRMSE_2D_radial(results_radial_ed, result_radial_ref, timesteps_radial)
        rad_ned_mean, rad_ned_max = self.__getRMSE_2D_radial(results_radial_ned, result_radial_ref, timesteps_radial)

        # Plotting
        self.__format_plot()

        fig, axs = plt.subplots(1, 2, figsize=(10.5, 5), constrained_layout=True)
        colors = self.__get_colors()

        n_axials = []
        # Plot axial equidistant
        for i, result in enumerate(results_axial_ed):
            n_axial = len(result.reactor.axial_discretization.get_centroids())
            n_axials.append(n_axial)
            axs[0].plot(n_axial, ax_ed_mean[i], marker='v', markersize=10, color=colors[0])
            axs[0].plot(n_axial, ax_ed_max[i], marker='x', markersize=10, color=colors[0])

        # Plot axial non-equidistant
        for i, result in enumerate(results_axial_ned):
            n_axial = len(result.reactor.axial_discretization.get_centroids())
            axs[0].plot(n_axial, ax_ned_mean[i], marker='v', markersize=10, color=colors[1])
            axs[0].plot(n_axial, ax_ned_max[i], marker='x', markersize=10, color=colors[1])

        max_n_rounded = math.ceil(np.max(n_axials) / 10) * 10
        axs[0].set_xlim(-max_n_rounded*0.05, max_n_rounded*1.05)
        axs[0].set_xticks(np.linspace(0, max_n_rounded, 5))
        axs[0].set_ylabel("normalized RMSE / -")
        axs[0].set_xlabel("number of axial volumes/ -")
        n_ref = len(result_axial_ref.reactor.axial_discretization.get_centroids())
        axs[0].set_title("reference: "+  str(n_ref) + " axial volumes" )

        n_radials = []
        # Plot radial equidistant
        for i, result in enumerate(results_radial_ed):
            n_radial = len(result.reactor.radial_discretization.get_centroids())
            n_radials.append(n_radial)
            axs[1].plot(n_radial, rad_ed_mean[i], marker='v', markersize=10, color=colors[0])
            axs[1].plot(n_radial, rad_ed_max[i], marker='x', markersize=10, color=colors[0])

        # Plot radial non-equidistant
        for i, result in enumerate(results_radial_ned):
            n_radial = len(result.reactor.radial_discretization.get_centroids())
            axs[1].plot(n_radial, rad_ned_mean[i], marker='v', markersize=10, color=colors[1])
            axs[1].plot(n_radial, rad_ned_max[i], marker='x', markersize=10, color=colors[1])

        max_n_rounded = math.ceil(np.max(n_radials) / 10) * 10
        axs[1].set_xlim(-max_n_rounded*0.05, max_n_rounded*1.05)
        axs[1].set_xticks(np.linspace(0, max_n_rounded, 5))
        axs[1].set_xlabel("number of radial volumes/ -")
        n_ref = len(result_radial_ref.reactor.radial_discretization.get_centroids())
        n_ref_ax = len(result_radial_ref.reactor.axial_discretization.get_centroids())
        axs[1].set_title("reference: " + str(n_ref_ax) + " axial, "+ str(n_ref) + " radial volumes")

        # set legend
        plt.plot([], [], marker='*', markersize=10, color="black", linestyle="None", label="maximum")
        plt.plot([], [], marker='v', markersize=10, color="black", linestyle="None",label="mean")
        plt.plot([], [], linestyle="None",label="equidistant")
        plt.plot([], [], linestyle="None",label="non-equidistant")
        legend = fig.legend(loc="upper center",
                   bbox_to_anchor=(0.5, 1.12), ncol=4)
        for text in legend.get_texts():
            if text.get_text() == 'equidistant':
                text.set_color(colors[0])
            elif text.get_text() == 'non-equidistant':
                text.set_color(colors[1])

        axs[0].yaxis.set_major_locator(MaxNLocator(nbins=3))  # Approx. ticks on y-axis
        axs[1].yaxis.set_major_locator(MaxNLocator(nbins=3))  # Approx. ticks on y-axis

        self.__export(foldername, "discretization_study", plt)
        plt.show()

    ####### additional methods for discretization study ########

    def __getRMSE_1D(self, results, result_ref, timesteps):
        start, stop = 0, results[0].reactor.reactorLength
        ## Calculating RMSEs for every Result
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
                    mean_nRMSEs_variables.append(self.__nRMSE(ref_values_interp, values_interp))

                mean_nRMSEs_timesteps.append(np.mean(mean_nRMSEs_variables))
                max_nRMSEs_timesteps.append(np.max(mean_nRMSEs_variables))
            mean_nRMSEs.append(np.mean(mean_nRMSEs_timesteps))
            max_nRMSEs.append(np.max(max_nRMSEs_timesteps))

        return mean_nRMSEs, max_nRMSEs

    def __getRMSE_2D_radial(self, results, result_ref, timesteps):
        start, stop = 0, results[0].reactor.reactorDiameter / 2
        ## Calculating RMSEs for every Result
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
                        mean_nRMSEs_variables.append(self.__nRMSE(ref_values_interp, values_interp))

                    mean_nRMSEs_axials.append(np.mean(mean_nRMSEs_variables))
                    max_nRMSEs_axials.append(np.max(np.percentile(mean_nRMSEs_variables, 90)))
                mean_nRMSEs_timesteps.append(np.mean(mean_nRMSEs_axials))
                max_nRMSEs_timesteps.append(np.mean(max_nRMSEs_axials))
            mean_nRMSEs.append(np.mean(mean_nRMSEs_timesteps))
            max_nRMSEs.append(np.mean(max_nRMSEs_timesteps))

        return mean_nRMSEs, max_nRMSEs

    def __nRMSE(self, values_ref, values):
        RMSE = self.__RMSE(values_ref, values)
        return RMSE / (np.mean(values))

    def __RMSE(self, values_ref, values):
        RMSE = np.sqrt(1 / len(values) * np.sum((values - values_ref) ** 2))
        return RMSE

    def __interpolate_arrays(self, values, centroids, start, stop):
        values = np.array(values)
        centroids = np.array(centroids)

        valid_mask = np.isfinite(centroids) & np.isfinite(values)
        centroids = centroids[valid_mask]
        values = values[valid_mask]

        resolution = 50
        common_orts = np.linspace(start, stop, resolution)

        interp1 = interp1d(centroids, values, kind='linear', fill_value='extrapolate')

        return np.column_stack((common_orts, interp1(common_orts)))

    '''################################ PLOT METHODS  - xxxxxxxxxxxxxxxxx ################################'''

    def plot_ignitionArc(self, foldername, results_ignition, results_extinction, T_wall_ign, T_wall_ext, timestep):
        self.__format_plot()
        dim = results_ignition[0].reactor.dimension

        fig, axs = plt.subplots(1, 1, figsize=(5.5, 5), constrained_layout=True)
        colors = self.__get_colors()

        # Plot Bremer Data
        T_bremer_ignition = np.array([311, 323, 342, 359, 377, 401, 424, 442, 454, 460, 465, 490, 512, 537, 561])
        X_bremer_ignition = np.array([0.00333, 0.00355, 0.00388, 0.00418, 0.00776, 0.00980, 0.0298, 0.0871, 0.201,
                             0.392, 0.930, 0.965, 0.980, 0.985, 0.979])

        T_bremer_extinction = np.array([308, 322, 336, 348, 359, 365, 371, 377, 388, 407, 425, 449, 471, 491, 513, 531, 549])
        X_bremer_extinction = np.array([0.00329, 0.00353, 0.00214, 0.00399, 0.00417, 0.686, 0.826, 0.845, 0.860,
                     0.882, 0.902, 0.931, 0.958, 0.971, 0.985, 0.985, 0.984])

        axs.plot(T_bremer_ignition, X_bremer_ignition, color="grey", linestyle="-", linewidth=2.5)
        axs.plot(T_bremer_extinction, X_bremer_extinction, color="grey", linestyle="-", linewidth=2.5)


        # Plot Simulation Data
        X_CO2_ign = []
        for result in results_ignition:
            if dim == 1:
                x_CO2_ign = result.getConversion_1D(timestep, 2)
            else:
                x_CO2_ign = result.getConversion_2D(timestep, 2)
            X_CO2_ign.append(result.average_trapezoidal(x_CO2_ign[-1, :]))

        X_CO2_ext = []
        for result in results_extinction:
            if dim == 1:
                x_CO2_ext = result.getConversion_1D(timestep, 2)
            else:
                x_CO2_ext = result.getConversion_2D(timestep, 2)
            X_CO2_ext.append(result.average_trapezoidal(x_CO2_ext[-1, :]))


        axs.plot(T_wall_ext, X_CO2_ext, color=colors[0], linestyle="-", linewidth=2.5, marker='v',
                 markersize=6, markerfacecolor="black", markeredgecolor='black', label = "extinction")
        axs.plot(T_wall_ign, X_CO2_ign, color=colors[1], linestyle="--", linewidth = 2.5,  marker='v',
                 markersize=6, markerfacecolor="black",markeredgecolor='black', label = "ignition")

        # Axis
        axs.set_xlabel(r'$T_{\mathregular{wall}} / K$')
        axs.set_ylabel(r'$X_{CO_2}/ -$')

        axs.set_xlim(275, 625)
        axs.set_xticks([300, 450, 600])

        axs.set_yticks([0, 0.5, 1])
        axs.xaxis.set_major_locator(MaxNLocator(nbins=3))  # Approx. ticks on y-axis

        #set legend for axs
        axs.legend(fontsize="16", loc="lower right")

        # # set title
        # axs.set_title("ignition and extinction arcs\n\n")

        # set legend for fig
        h1, = plt.plot([], [], linestyle = "-", markersize=10, color="grey", label="Bremer et. al.", linewidth=2.5)
        h2, = plt.plot([], [], marker='v', markersize=10, color="black", linestyle="-", label="simulation", linewidth=2.5)
        fig.legend(handles=[h1, h2], loc="upper center", bbox_to_anchor=(0.58, 1.1), ncol=2, fontsize="16")

        self.__export(foldername, "arcs_2D", plt)
        plt.show()






































    def plot1D(self, name, result, timestep):

        w_i, T, u, p = result.get_values(timestep)
        z_pos = result.get_z_pos() / result.reactor.reactorLength

        w_f = CasADi.SX.sym('w_f', 4)
        T_f = CasADi.SX.sym('T_f')
        p_f = CasADi.SX.sym('p_f')
        u_f = CasADi.SX.sym('u_f')

        Reynolds = CasADi.Function('f_Re_CasADi', [w_f, T_f, u_f, p_f], [result.reactor.Re(w_f, T_f,u_f, p_f)])

        check = Reynolds(w_i[z_pos, timestep, :], T[z_pos, timestep], p[z_pos, timestep], u[z_pos, timestep])

        fig, ax = plt.subplots(1, 1, figsize=(4.2, 5.7), constrained_layout=True, sharex=True)

        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        ax.plot(z_pos, check, color=colors[0], linestyle="--")



    def plot_Twall_vs_validation(self, name, result, timestep):
        fig, axs = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        # Plot Simulation Data Wall & Center
        w_i, T, p, u = result.get_2D_values(timestep)
        z_pos = result.get_z_pos() / result.reactor.reactorLength

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





    def plot1D_vsPseudo1D_vs_val(self, name, result1D, resultP1D, timestep):
        fig, axs = plt.subplots(4, 1, figsize=(4.2, 5.7), constrained_layout=True, sharex=True)
        colors = plt.cm.Dark2(np.linspace(0, 1, 8))

        # Plot Simulation Data 1D
        w_i, T, p, u = result1D.get_values(timestep)
        z_pos = result1D.get_z_pos() / result1D.reactor.reactorLength

        axs[0].plot(z_pos, w_i[:, 0], color=colors[0], linestyle="--", label="CH4")
        axs[0].plot(z_pos, w_i[:, 1], color=colors[1], linestyle="--", label="H2O")
        axs[0].plot(z_pos, w_i[:, 2], color=colors[2], linestyle="--", label="CO2")
        axs[0].plot(z_pos, w_i[:, 3], color=colors[3], linestyle="--", label="H2")
        axs[1].plot(z_pos, T, color=colors[0], linestyle="--")
        axs[2].plot(z_pos, u, color=colors[0], linestyle="--")
        axs[3].plot(z_pos, p * 1e-5, color=colors[0], linestyle="--", label="1D")

        # Plot Simulation Data P1D
        w_i_2D, T_2D, p_2D, u_2D = resultP1D.get_values(timestep)
        z_pos_2D = resultP1D.get_z_pos() / resultP1D.reactor.reactorLength

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


    def plot2D_T_half(self, name, result, timestep):
        self.__format_plot()
        fig, axs = plt.subplots(1, 1, figsize=(20, 5))

        T = result.get_2D_values(timestep)[1]

        z_pos = result.get_z_pos() / result.reactor.reactorLength
        r_pos = 2 * result.get_r_pos() / result.reactor.reactorDiameter

        tol_cmap = TOLcmaps()
        cmp = tol_cmap.get('sunset')

        z_mesh, r_mesh = np.meshgrid(z_pos, r_pos, indexing='ij')

        T_min, T_max = self.__get_min_max(T)

        axs.pcolormesh(z_mesh, r_mesh, T, cmap=cmp, vmin=T_min, vmax=T_max, shading='gouraud')

        axs.set_axis_off()  # Hides everything: axes, ticks, and labels
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)  # Removes padding

        plt.savefig("results/T_half_plot.svg", format="svg", bbox_inches="tight")
        plt.savefig("results/T_half_plot.pdf", format="pdf", bbox_inches="tight")
        plt.savefig("results/T_half_plot.eps", format="eps", bbox_inches="tight")
        plt.savefig("results/T_half_plot.png", format="png", dpi=600, bbox_inches="tight")

        plt.show()
