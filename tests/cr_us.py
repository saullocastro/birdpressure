from BirdPressure.Bird import Bird
from BirdPressure.ImpactScenario import ImpactScenario
from BirdPressure.EOS import EOS
from BirdPressure.Timing import Timing

import matplotlib.pyplot as plt
import numpy as np


bird = Bird()
scen_range = ImpactScenario(np.linspace(0.1,300.1,31),90,False)
timing_range = Timing(bird,scen_range)
PH_range = timing_range.det_peak_P()

EOS_lin = EOS(bird, 'Linear_EOS')
EOS_mur = EOS(bird, 'Murnaghan')
EOS_pol = EOS(bird, 'Polynomial_EOS')
lin_cr_range = np.zeros_like(PH_range)
mur_cr_range = np.zeros_like(PH_range)
pol_cr_range = np.zeros_like(PH_range)

for idx,PH in enumerate(PH_range):
    lin_cr_range[idx] = EOS_lin.get_c_r(PH)
    mur_cr_range[idx] = EOS_mur.get_c_r(PH)
    pol_cr_range[idx] = EOS_pol.get_c_r(PH)


cm = 1/2.54

fig,ax = plt.subplots(figsize = (30*cm,20*cm))
ax.plot(scen_range.get_impact_velocity(),timing_range.u_s, label = r'$u_s$')
ax.plot(scen_range.get_impact_velocity(),lin_cr_range, label = r'Linear EOS $c_r$')
ax.plot(scen_range.get_impact_velocity(),mur_cr_range, label = r'Murnaghan EOS $c_r$')
ax.plot(scen_range.get_impact_velocity(),pol_cr_range, label = r'Polynomial EOS $c_r$')
ax.legend(fontsize= 'xx-large', loc = 'upper left')
ax.grid()
ax.set_xlabel(r'Impact velocity, $u_0$ [m/s]',fontsize= 'xx-large')
ax.set_ylabel('Wave speed,[m/s]',fontsize= 'xx-large')
ax.set_xlim(xmin = 0)
ax.set_ylim(ymin = 0)
plt.show()