from BirdPressure.Bird import Bird
from BirdPressure.ImpactScenario import ImpactScenario

import matplotlib.pyplot as plt
import numpy as np

class Timing():
    def __init__(self, bird: Bird, scenario: ImpactScenario):
        self.bird = bird
        self.scenario = scenario


        self.u_s = self.find_u_s()
        self.t_D = self.find_t_D()


    def find_u_s(self):
        u_s = self.bird.get_sound_speed() + self.bird.get_k() * self.scenario.get_impact_velocity()
        return u_s

    def find_t_D(self):
        return self.bird.get_length()/self.scenario.get_impact_velocity()
    def find_t_c(self, release_wave_velocity):

        if self.check_release_wave_velocity(release_wave_velocity):
            return self.bird.get_length()/(np.sqrt(release_wave_velocity**2-(self.u_s-self.scenario.get_impact_velocity())**2))
    def find_t_b(self,release_wave_velocity):
        return self.bird.get_diameter()/(2 * release_wave_velocity )

    def critical_aspect_ratio(self,release_wave_velocity):

        if self.check_release_wave_velocity(release_wave_velocity):
            return self.u_s / (2* np.sqrt(release_wave_velocity ** 2 - (self.u_s - self.scenario.get_impact_velocity()) ** 2))

    def det_peak_P(self, use_us = True, use_v_norm = True):
        if use_v_norm:
            velocity = self.scenario.get_normal_velocity()
        else:
            velocity = self.scenario.get_impact_velocity()
        if use_us:
            return self.bird.get_density()*self.u_s*velocity
        return self.bird.get_density()*self.bird.get_sound_speed()*velocity

    def check_release_wave_velocity(self,release_wave_velocity):
        if type(release_wave_velocity) == int or type(release_wave_velocity) == float or type(release_wave_velocity) == np.float64:

            if (self.u_s - self.scenario.get_impact_velocity()) > release_wave_velocity:
                print('Error, u_s is much higher then c_r. Use different EOS')
                return False

            else:
                return True

        else:
            for idx, cr in enumerate(release_wave_velocity):
                if (self.u_s - self.scenario.get_impact_velocity())[idx] > cr:
                    print('Error, u_s is much higher then c_r. Use different EOS')
                    return False
                    break
                else:
                    return True


    def find_force_history(self, f_peak, f_steady, release_wave_velocity):
        self.t_b = self.find_t_b(release_wave_velocity)

        f_history = np.array([f_peak, f_peak,f_steady,f_steady])
        t_history = np.array([0,self.t_b,self.t_b, self.t_D])

        j_tot = f_peak * self.t_b + f_steady * (self.t_D-self.t_b)

        return f_history, t_history, j_tot

    def plot_f_history(self,f_hist,t_hist):
        cm = 1 / 2.54
        fig, ax = plt.subplots(figsize=(30 * cm, 20 * cm))

        ax.plot(t_hist* 10 ** (3), f_hist * 10 ** (-3))
        ax.grid()
        ax.set_xlabel(r'Time, [\mu s]', fontsize='xx-large')
        ax.set_ylabel(r'Average Force [kN]', fontsize='xx-large')
        ax.set_xlim(xmin=0)
        ax.set_ylim(ymin=0)
        plt.show()



