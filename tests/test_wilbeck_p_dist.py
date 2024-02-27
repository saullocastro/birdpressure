from birdpressure.Bird import Bird
from birdpressure.ImpactScenario import ImpactScenario
from birdpressure.GridGenerator import GridGenerator
from birdpressure.EOS import EOS
from birdpressure.AnalWilbeck import AnalWilbeck

import matplotlib.pyplot as plt

def test_wilbeck_p_dist():
    bird_1 = Bird()
    scenario = ImpactScenario(impact_velocity= 100, impact_angle=20, use_radian = False )
    grid = GridGenerator(bird_1,scenario,21)
    eos = EOS(bird_1, 'Linear_EOS')


    wily = AnalWilbeck(bird_1,scenario, grid, eos)
    wily.show_banks()
    wily.show_leach()

    cm = 1 / 2.54

    fig, ax = plt.subplots(figsize=(30 * cm, 20 * cm))
    ax.plot(wily.r_leach/(bird_1.get_diameter()/2), wily.P_leach * 10 ** (-6), label='Leach and Walker')
    ax.plot(wily.r_banks/(bird_1.get_diameter()/2), wily.P_banks * 10 ** (-6), label='Banks and Chandrasekhara')
    ax.grid()
    ax.set_xlabel(r'Distance from Impact center, $r/a$ [-]', fontsize='xx-large')
    ax.set_ylabel(r'Pressure, $P$ [MPa]', fontsize='xx-large')
    ax.set_xlim(xmin = 0, xmax = 3.5)
    ax.set_ylim(ymin = 0)
    ax.legend()
    plt.show()

