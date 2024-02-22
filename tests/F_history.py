from BirdPressure.Bird import Bird
from BirdPressure.ImpactScenario import ImpactScenario
from BirdPressure.GridGenerator import GridGenerator
from BirdPressure.EOS import EOS
from BirdPressure.AnalBarber import AnalBarber
from BirdPressure.AnalWilbeck import AnalWilbeck



bird_1 = Bird()
scenario = ImpactScenario(impact_velocity= 100, impact_angle=60, use_radian = False )
grid = GridGenerator(bird_1,scenario,31)
eos = EOS(bird_1, 'Polynomial_EOS')

barber = AnalBarber(bird_1,scenario, grid, eos)
wily = AnalWilbeck(bird_1, scenario, grid, eos)

barber.show_averaged_force_history()
wily.show_averaged_force_history()