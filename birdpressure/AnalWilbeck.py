from BirdPressure.Bird import Bird
from BirdPressure.GridGenerator import GridGenerator
from BirdPressure.ImpactScenario import ImpactScenario
from BirdPressure.Timing import Timing
from BirdPressure.EOS import EOS

import matplotlib.pyplot as plt
import numpy as np

class AnalWilbeck(Timing):

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

    # def det_peak_P(self, use_us = True):
    #     if use_us:
    #         return self.bird.get_density()*self.u_s*self.impact_scenario.get_impact_velocity()
    #     return self.bird.get_density()*self.bird.get_sound_speed()*self.impact_scenario.get_impact_velocity()

    def get_averages(self, Circular_cross = True):
        if Circular_cross:
            impact_A = np.pi * (self.bird.diameter/2)**2
        else:
            impact_A = np.pi * self.grid.min_axis * self.grid.maj_axis

        P_avg = self.bird.get_density() * self.impact_scenario.get_impact_velocity()**2 * np.sin(self.impact_scenario.get_impact_angle())
        F_avg = P_avg * impact_A

        return P_avg, F_avg

    def get_P_s(self):
        return 1/2 * self.bird.get_density() * self.impact_scenario.get_impact_velocity()**2 * np.sin(self.impact_scenario.get_impact_angle())
    def find_leach(self,r,xi_2 =2.58 ):
        P = self.P_s * (1 - 3 * (r / (xi_2 * (self.bird.get_diameter()/2))) ** 2 + 2 * (r / (xi_2 * (self.bird.get_diameter()/2))) ** 3)
        return P

    def find_banks(self,r,xi_1 = 0.5):
        P = self.P_s * np.e ** (-xi_1 * (r / (self.bird.get_diameter()/2)) ** 2)
        return P


    def show_leach(self):
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
        self.p_steady_av, self.f_steady_av = self.get_averages()
        self.f_peak_av = self.P_H*self.grid.impact_area
        self.f_hist, self.t_hist, self.impulse_tot = self.find_force_history(self.f_peak_av, self.f_steady_av, self.c_r)
        self.plot_f_history(self.f_hist, self.t_hist)



















