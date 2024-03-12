from birdpressure.Timing import Timing
from birdpressure.Bird import Bird
from birdpressure.GridGenerator import GridGenerator
from birdpressure.ImpactScenario import ImpactScenario
from birdpressure.EOS import EOS

import matplotlib.pyplot as plt
import numpy as np
class AnalBarber(Timing):
    """ Class for Wilbeck's analytical bird strike solution
        Parameters
        ----------
        bird: class
        impact_scenario: class
        grid: class
        eos: class
        initial_pressure: int, default: 101325
            ambient pressure default at sea-level

        Attributes
        ----------
        bird: class
        impact_scenario: class
        grid: class
        eos: class
        use_us: bool, default = True
            use shock velocity to determine peak pressure
        use_v_norm: bool, default = True
            use normal velocity component to determine peak pressure
        initial_pressure: int, default: 101325
            ambient pressure default at sea-level
        U: ndarray
            superimposed velocity in x direction [m/s]
        V: ndarray
            superimposed velocity in y direction [m/s]
        W: ndarray
            superimposed velocity in z direction [m/s]
        c_p: ndarray
            array with pressure coefficient over impact area [-]
        P: ndarray
            array with pressures over impact area [Pa]
        f_steady_av: float or ndarray
            average force exerted during steady state phase [N]
        P_steady_av: float or ndarray
            average force exerted during steady state phase [Pa]
        P_H: float or ndarray
            peak pressure (hugonoit pressure) [Pa]
        f_peak_av: float or ndarray
            force applied during hugonoit state [N]
        c_r: float or ndarray
            release wave velocity [m/s]
        t_b: float or ndarray
            peak impact pressure duration [s]
        t_c: float or ndarray
            duration until shocked area dissipates [s]
        f_hist: ndarray
            pressure history force array [N]
        t_hist: ndarray
            pressure history time array [s]
        impulse_tot: float or ndarray
            total impulse over per impact [Ns]
    """

    def __init__(self, bird: Bird, impact_scenario: ImpactScenario, grid: GridGenerator, eos: EOS, initial_pressure = 101325, use_us = True, use_v_norm = True):
        Timing.__init__(self, bird, impact_scenario, use_us, use_v_norm)

        
        self.initial_pressure = initial_pressure
        self.bird = bird
        self.impact_scenario = impact_scenario
        self.grid = grid
        self.eos = eos
        self.use_us = use_us
        self.use_v_norm = use_v_norm

        self.U, self.V, self.W = self.get_UVW()
        self.c_p = 1 - (self.V ** 2 + self.W ** 2) / (self.impact_scenario.get_impact_velocity() ** 2)
        self.P = self.c_p * 1 / 2 * self.bird.get_density() * self.impact_scenario.get_impact_velocity() ** 2
        self.f_steady_av, self.P_steady_av = self.get_averages()
        self.P_H = self.det_peak_P()
        self.f_peak_av = self.P_H * np.pi * (self.bird.diameter/2)**2
        self.c_r = self.eos.get_c_r(self.P_H)
        self.t_b = self.find_t_b(self.c_r)
        self.t_c = self.find_t_c(self.c_r)
        self.f_hist, self.t_hist, self.impulse_tot = self.find_force_history(self.f_peak_av, self.f_steady_av, self.c_r)



    def get_UVW(self):
        """Function to compute superimposd velocity in x,y and z direction."""
        uvw = np.zeros([self.grid.ids.shape[0], 3])
        for idx_1, id_cent in enumerate(self.grid.ids):

            x = 0
            y = self.grid.centers[id_cent[0], 0]
            z = self.grid.centers[id_cent[1], 1]

            for idx in np.arange(self.grid.corners.shape[0]):
                eta_1 = self.grid.corners[idx, 0]
                eta_2 = self.grid.corners[idx, 1]
                zeta_1 = self.grid.corners[idx, 2]
                zeta_2 = self.grid.corners[idx, 3]

                u, v, w, r_1, r_2, r_3, r_4 = self.get_velocity(x, y, z, eta_1, eta_2, zeta_1, zeta_2)
                uvw[idx_1, :] += [u, v, w]

        U = self.impact_scenario.get_impact_velocity() * np.sin(self.impact_scenario.get_impact_angle()) + \
            self.impact_scenario.get_impact_velocity() * np.sin(self.impact_scenario.get_impact_angle()) / (2 * np.pi) * uvw[:, 0]

        V = self.impact_scenario.get_impact_velocity() * np.cos(self.impact_scenario.get_impact_angle()) + \
            self.impact_scenario.get_impact_velocity() * np.sin(self.impact_scenario.get_impact_angle()) / (2 * np.pi) * uvw[:, 1]

        W = self.impact_scenario.get_impact_velocity() * np.sin(self.impact_scenario.get_impact_angle()) / (2 * np.pi) * uvw[:, 2]

        return U, V, W


    def get_velocity(self, x, y, z, eta_1, eta_2, zeta_1, zeta_2):
        """Function to compute velocity field induced by uniformly distributed sources."""

        r_1 = np.sqrt(x ** 2 + (y - eta_1) ** 2 + (z - zeta_1) ** 2)
        r_2 = np.sqrt(x ** 2 + (y - eta_2) ** 2 + (z - zeta_1) ** 2)
        r_3 = np.sqrt(x ** 2 + (y - eta_2) ** 2 + (z - zeta_2) ** 2)
        r_4 = np.sqrt(x ** 2 + (y - eta_1) ** 2 + (z - zeta_2) ** 2)

        if x == 0:
            u = (np.pi / 2 * (((z - zeta_2) * (y - eta_2)) / abs((z - zeta_2) * (y - eta_2))) + \
                 np.pi / 2 * ((z - zeta_1) * (y - eta_1)) / abs((z - zeta_1) * (y - eta_1)) - \
                 np.pi / 2 * ((z - zeta_1) * (y - eta_2) / abs((z - zeta_1) * (y - eta_2))) - \
                 np.pi / 2 * ((z - zeta_2) * (y - eta_1) / abs((z - zeta_2) * (y - eta_1))))
        else:

            u = (np.arctan((z - zeta_2) * (y - eta_2) / (x * r_3)) + np.arctan((z - zeta_1) * (y - eta_1) / (x * r_1)) - \
                 np.arctan((z - zeta_1) * (y - eta_2) / (x * r_2)) - np.arctan((z - zeta_2) * (y - eta_1) / (x * r_4)))

        v = np.log(((r_3 + (zeta_2 - z)) * (r_1 + (zeta_1 - z))) / ((r_4 + (zeta_2 - z)) * (r_2 + (zeta_1 - z))))
        w = np.log(((r_3 + (eta_2 - y)) * (r_1 + (eta_1 - y))) / ((r_2 + (eta_2 - y)) * (r_4 + (eta_1 - y))))

        return u, v, w, r_1, r_2, r_3, r_4

    def mesh_results(self, result):
        """Function to save results in mesh format."""

        result_cells = np.zeros_like(self.grid.mesh_X)

        for idx, id in enumerate(self.grid.ids):
            result_cells[id[0], id[1]] = result[idx]

        return result_cells

    def plot_result(self, results, results_name: str, result_unit: str):
        """Function to plot results"""

        result_cells = self.mesh_results(results)

        cm = 1 / 2.54

        fig, ax = plt.subplots(figsize=(20 * cm, 20 * cm))

        ax.set_aspect('equal')
        ax.set_title(f""""Contour Plot of {results_name}""")
        cf = ax.pcolormesh(self.grid.mesh_X, self.grid.mesh_Y,result_cells)
        fig.colorbar(cf, ax=ax, label=result_unit)

        plt.show()

    def get_axis_data(self,result):
        """Function to plot results along major and minor axis

        Notes
        -----
        In order to use this function a grid with uneven number of elements needs to be used.
        """
        if (self.grid.n_elements%2) == 0:
            print('Error: n_elements needs to be odd to generate axis data')
        else:
            id_major = []
            id_minor = []
            maj_ax_val = []
            min_ax_val = []

            for id in self.grid.ids:
                y = round(self.grid.centers[id[0], 0], 9)
                z = round(self.grid.centers[id[1], 1], 9)

                if z == 0:
                    id_major.append(id)
                    maj_ax_val.append(y)
                if y == 0:
                    id_minor.append(id)
                    min_ax_val.append(z)

            id_major = np.array(id_major)
            id_minor = np.array(id_minor)
            maj_ax_val = np.array(maj_ax_val)
            min_ax_val = np.array(min_ax_val)

            maj_result = np.zeros_like(maj_ax_val)
            min_result = np.zeros_like(min_ax_val)

            result_cells = self.mesh_results(result)

            for n,id in enumerate(id_major):
                maj_result[n] = result_cells[id[0],id[1]]

            for n,id in enumerate(id_minor):
                min_result[n] = result_cells[id[0],id[1]]

            return maj_result, maj_ax_val, min_result, min_ax_val

    def get_averages(self):
        """Function to compute steady state average values."""

        cell_area = self.grid.resolution**2
        impact_area = cell_area*len(self.P)
        force = 0

        for pressure in self.P:
            force+= pressure*cell_area

        average_P = force/impact_area

        return force, average_P

    def show_averaged_force_history(self):
        """Function to plot force history"""
        self.plot_f_history(self.f_hist,self.t_hist)







