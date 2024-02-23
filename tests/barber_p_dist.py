from birdpressure.Bird import Bird
from birdpressure.ImpactScenario import ImpactScenario
from birdpressure.GridGenerator import GridGenerator
from birdpressure.EOS import EOS
from birdpressure.AnalBarber import AnalBarber



bird_1 = Bird()
scenario = ImpactScenario(impact_velocity= 100, impact_angle=90, use_radian = False )
grid = GridGenerator(bird_1,scenario,21)
eos = EOS(bird_1, 'Polynomial_EOS')

barber = AnalBarber(bird_1,scenario, grid, eos)

barber.plot_result(barber.P *10**(-6), 'Steady State Pressure', '[MPa]')
barber.plot_result(barber.c_p,'Steady State a Pressure Coefficient', '[-]')
