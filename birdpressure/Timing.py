from birdpressure.Bird import Bird
from birdpressure.ImpactScenario import ImpactScenario

import matplotlib.pyplot as plt
import numpy as np

class Timing():
    """Compute impact timing as well as force history

    Notes
    -----
    This class is inherited by AnalWilbeck and AnalBarber to provide timing functions/methods

    Parameters
    ----------
    bird: class
    impact_scenario: class

    Attributes
    ----------
    bird: class
    impact_scenario: class
    u_s: float or ndarray
        shock wave velocity [m/s]
    t_D: float or ndarray
        total impact duration (squash time) [s]
    """

    def __init__(self, bird: Bird, impact_scenario: ImpactScenario):
        self.bird = bird
        self.impact_scenario = impact_scenario


        self.u_s = self.find_u_s()
        self.t_D = self.find_t_D()


    def find_u_s(self):
        """Linear hugonoit function to find shock velocity.

        Returns
        -------
        u_s: float or ndarray
            shock wave velocity [m/s]
        """
        u_s = self.bird.get_sound_speed() + self.bird.get_k() * self.impact_scenario.get_impact_velocity()
        return u_s

    def find_t_D(self):
        """Compute total impact duration.

        Returns
        -------
        t_D: float or ndarray
            total impact duration (squash time) [s]
        """
        return self.bird.get_length()/self.impact_scenario.get_impact_velocity()

    def find_t_c(self, release_wave_velocity):
        """Compute shock state duration.
            Parameters
            ----------
            release_wave_velocity: float or ndarray
                release wave velocity (c_r) P vs rho slope at hugonoit pressure [m/s]

            Returns
            -------
            t_c: float or ndarray
                duration until shocked area dissipates [s]
        """

        if self.check_release_wave_velocity(release_wave_velocity):
            self.t_c = self.bird.get_length()/(np.sqrt(release_wave_velocity**2-(self.u_s-self.impact_scenario.get_impact_velocity())**2))
            return self.t_c
    def find_t_b(self,release_wave_velocity):
        """Compute peak impact pressure duration.

            Parameters
            ----------
            release_wave_velocity: float or ndarray
                release wave velocity (c_r) P vs rho slope at hugonoit pressure [m/s]

            Returns
            -------
            t_b: float or ndarray
                peak impact pressure duration [s]
            """
        self.t_b = self.bird.get_diameter()/(2 * release_wave_velocity )
        return self.t_b

    def critical_aspect_ratio(self,release_wave_velocity):
        """ Compute critical aspect ratio

            Notes
            -----
            May return nothing if release wave velocity is too low

            Parameters
            ----------
            release_wave_velocity: float or ndarray
                release wave velocity (c_r) P vs rho slope at hugonoit pressure [m/s]

            Returns
            -------
            critical_aspect: float or ndarray
                critical aspect ratio (L/D) [-]
        """

        if self.check_release_wave_velocity(release_wave_velocity):
            critical_aspect = self.u_s / (2* np.sqrt(release_wave_velocity ** 2 - (self.u_s - self.impact_scenario.get_impact_velocity()) ** 2))
            return critical_aspect

    def det_peak_P(self, use_us = True, use_v_norm = True):
        """ Compute peak pressure (hugonoit pressure)

            Parameters
            ----------
            use_us: bool, default: True
                decide whether to use shock velocity or speed of sound
            use_v_norm: bool, default: True
                decide whether to use normal velocity or impact velocity

            Returns
            -------
            P_H: float or ndarray
                peak pressure [N/m^2]

        """
        if use_v_norm:
            velocity = self.impact_scenario.get_normal_velocity()
        else:
            velocity = self.impact_scenario.get_impact_velocity()
        if use_us:
            P_H = self.bird.get_density()*self.u_s*velocity
        else:
            P_H = self.bird.get_density()*self.bird.get_sound_speed()*velocity
        return P_H
    def check_release_wave_velocity(self,release_wave_velocity):
        """Function to check whether release velocity (c_r) is high enough

            Parameters
            ----------
            release_wave_velocity: float or ndarray
                release wave velocity (c_r) P vs rho slope at hugonoit pressure [m/s]

            Returns
            -------
            bool
                True if high enough, False if not
        """
        check = (self.u_s - self.impact_scenario.get_impact_velocity())

        if isinstance(check,np.ndarray) or isinstance(release_wave_velocity,np.ndarray):
            if (check > release_wave_velocity).all():
                print('Error, u_s is much higher then c_r. Use different EOS.')
                return False
            else:
                return True

        else:
            if check > release_wave_velocity:
                print('Error, u_s is much higher then c_r. Use different EOS.')
                return False
            else:
                return True


    def find_force_history(self, f_peak, f_steady, release_wave_velocity):
        """ Create numpy arrays for force history

            Parameters
            ----------
            f_peak: float
                peak force [N]
            f_steady: float
                steady state force [N]
            release_wave_velocity: float
                release wave velocity (c_r) P vs rho slope at hugonoit pressure [m/s]

            Returns
            -------
            f_history: ndarray
                numpy array of force history forces
            t_history: ndarray
                numpy array of force history time steps
            j_tot: float
                total impulse [Ns]
        """
        self.t_b = self.find_t_b(release_wave_velocity)

        f_history = np.array([f_peak, f_peak,f_steady,f_steady,0])
        t_history = np.array([0,self.t_b,self.t_b, self.t_D, self.t_D])

        impulse_tot = f_peak * self.t_b + f_steady * (self.t_D-self.t_b)

        return f_history, t_history, impulse_tot

    def plot_f_history(self,f_hist,t_hist):
        """Function that plots force history """
        cm = 1 / 2.54
        fig, ax = plt.subplots(figsize=(30 * cm, 20 * cm))

        ax.plot(t_hist* 10 ** (3), f_hist * 10 ** (-3))
        ax.grid()
        ax.set_xlabel(r'Time, [$\mu s$]', fontsize='xx-large')
        ax.set_ylabel(r'Average Force [kN]', fontsize='xx-large')
        ax.set_xlim(xmin=0)
        ax.set_ylim(ymin=0)
        plt.show()



