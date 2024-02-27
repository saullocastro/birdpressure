from birdpressure.Bird import Bird
from birdpressure.ImpactScenario import ImpactScenario
from birdpressure.GridGenerator import GridGenerator
from birdpressure.EOS import EOS
from birdpressure.AnalBarber import AnalBarber

import numpy as np
import matplotlib.pyplot as plt

def test_barber_axial_p_dist():
    bird_1 = Bird()
    scenario = ImpactScenario(impact_velocity= 100, impact_angle=20, use_radian = False )
    grid = GridGenerator(bird_1,scenario,21)
    eos = EOS(bird_1, 'Polynomial_EOS')

    barber = AnalBarber(bird_1,scenario, grid, eos)

    cp_maj, maj_ax, cp_min, min_ax = barber.get_axis_data(barber.c_p)

    cm = 1 / 2.54
    fig,ax = plt.subplots(figsize = (20*cm,20*cm))
    ax.plot(maj_ax/np.max(maj_ax), cp_maj, color = 'green', label = 'Major Axis')
    ax.plot(min_ax/np.max(min_ax), cp_min, color = 'blue', label = 'Minor Axis')

    #ax.set_xlim(a, -a)
    #ax.set_ylim(ymin = -a, ymax = a)
    ax.set_xlabel('Nondimensional radius [-]')
    ax.set_ylabel('C_p [-]')
    ax.legend()
    plt.show()