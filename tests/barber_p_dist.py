from BirdPressure.Bird import Bird
from BirdPressure.ImpactScenario import ImpactScenario
from BirdPressure.GridGenerator import GridGenerator
from BirdPressure.EOS import EOS
from BirdPressure.AnalBarber import AnalBarber



bird_1 = Bird()
scenario = ImpactScenario(impact_velocity= 100, impact_angle=90, use_radian = False )
grid = GridGenerator(bird_1,scenario,21)
eos = EOS(bird_1, 'Polynomial_EOS')

barber = AnalBarber(bird_1,scenario, grid, eos)

barber.plot_result(barber.P *10**(-6), 'Steady State Pressure', '[MPa]')
barber.plot_result(barber.c_p,'Steady State a Pressure Coefficient', '[-]')
