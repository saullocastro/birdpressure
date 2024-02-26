from birdpressure.Bird import Bird
from birdpressure.GridGenerator import GridGenerator
from birdpressure.ImpactScenario import ImpactScenario
from birdpressure.Timing import Timing
from birdpressure.EOS import EOS

import matplotlib.pyplot as plt
import numpy as np

class AnalWilbeck(Timing):
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
        initial_pressure: int, default: 101325
            ambient pressure default at sea-level
        eos: class
        bird: class
        impact_scenario: class
        grid: class
        u_s: float or ndarray
            shock wave velocity [m/s]
        P_H: float or ndarray
            peak (hugonoit) pressure [N/m^2]
        c_r: float or ndarray
            release wave velocity [m/s]
        P_s: float ndarray
            steady state stagnation point pressure [N/m^2]
    """

    def __init__(self, bird: Bird, impact_scenario: ImpactScenario, grid: GridGenerator, eos: EOS,  initial_pressure = 101325):
        Timing.__init__(self, bird, impact_scenario)

        self.initial_pressure = initial_pressure
        self.eos = eos
        self.bird = bird
        self.impact_scenario = impact_scenario
        self.grid = grid

        self.u_s = self.find_u_s()
        self.P_H = self.det_peak_P()
        self.c_r = self.eos.get_c_r(self.P_H)
        self.P_s = self.get_P_s()

    def get_averages(self, circular_cross = True):
        """Function to compute average pressure and force during steady state"""

        if circular_cross:
            impact_A = np.pi * (self.bird.diameter/2)**2
        else:
            impact_A = np.pi * self.grid.min_axis * self.grid.maj_axis

        P_avg = self.bird.get_density() * self.impact_scenario.get_impact_velocity()**2 * np.sin(self.impact_scenario.get_impact_angle())
        F_avg = P_avg * impact_A

        return P_avg, F_avg

    def get_P_s(self):
        """Funtion to compute steady state stagnation point pressure"""
        return 1/2 * self.bird.get_density() * self.impact_scenario.get_impact_velocity()**2 * np.sin(self.impact_scenario.get_impact_angle())
    def find_leach(self,r,xi_2 =2.58 ):
        """Function for Leach & Walker steady state pressure distribution.
        Parameters
        ----------
        r: float or ndarray
            distance from center [m]
        xi_2: float, default: 2.58
            material constant [-]

        Returns
        -------
        P: float or ndarray
            Pressure distribution [Pa]
        """
        P = self.P_s * (1 - 3 * (r / (xi_2 * (self.bird.get_diameter()/2))) ** 2 + 2 * (r / (xi_2 * (self.bird.get_diameter()/2))) ** 3)
        return P

    def find_banks(self,r,xi_1 = 0.5):
        """Function for Banks & Chandrasekhara steady state pressure distribution.
                Parameters
                ----------
                r: float or ndarray
                    distance from center [m]
                xi_1: float, default: 0.5
                    material constant [-]

                Returns
                -------
                P: float or ndarray
                    Pressure distribution [Pa]
                """
        P = self.P_s * np.e ** (-xi_1 * (r / (self.bird.get_diameter()/2)) ** 2)
        return P


    def show_leach(self):
        """Funtion to plot Leach & Walker steady state pressure distribution """
        r = np.linspace(0,4*self.bird.get_length()/2,50)
        P_dist = self.find_leach(r)

        min_idx = np.argmin(P_dist)

        P_dist = P_dist[:min_idx+1]
        r      = r[:min_idx + 1]

        cm = 1 / 2.54

        fig, ax = plt.subplots(figsize=(30 * cm, 20 * cm))
        ax.plot(r, P_dist * 10 ** (-6), label='Leach and Walker')
        ax.grid()
        ax.set_xlabel(r'Distance from Impact center, $r$ [m]', fontsize='xx-large')
        ax.set_ylabel(r'Pressure, $P$ [MPa]', fontsize='xx-large')
        plt.show()

        self.r_leach = r
        self.P_leach = P_dist

    def show_banks(self):
        """Funtion to plot Banks & Chandrasekhara steady state pressure distribution """
        r = np.linspace(0,4*self.bird.get_length()/2,50)
        P_dist = self.find_banks(r)

        min_idx = np.argmin(P_dist)

        P_dist = P_dist[:min_idx + 1]
        r = r[:min_idx + 1]

        cm = 1 / 2.54

        fig, ax = plt.subplots(figsize=(30 * cm, 20 * cm))
        ax.plot(r, P_dist * 10 ** (-6), label='Banks and Chandrasekhara')
        ax.grid()
        ax.set_xlabel(r'Distance from Impact center, $r$ [m]', fontsize='xx-large')
        ax.set_ylabel(r'Pressure, $P$ [MPa]', fontsize='xx-large')
        plt.show()

        self.r_banks = r
        self.P_banks = P_dist


    def show_averaged_force_history(self):
        """Function to plot force history"""
        self.p_steady_av, self.f_steady_av = self.get_averages()
        self.f_peak_av = self.P_H*self.grid.impact_area
        self.f_hist, self.t_hist, self.impulse_tot = self.find_force_history(self.f_peak_av, self.f_steady_av, self.c_r)
        self.plot_f_history(self.f_hist, self.t_hist)



















