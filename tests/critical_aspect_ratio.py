from BirdPressure.Bird import Bird
from BirdPressure.ImpactScenario import ImpactScenario
from BirdPressure.EOS import EOS
from BirdPressure.Timing import Timing

import matplotlib.pyplot as plt
import numpy as np

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

u0_range = np.linspace(0.1,350.1,61)
bird = Bird()
scen_range = ImpactScenario(u0_range,90,False)
timing_range = Timing(bird,scen_range)
eos = EOS(bird, 'Polynomial_EOS')
PH_range = timing_range.det_peak_P()

cr_range = np.zeros_like(PH_range)

for idx,PH in enumerate(PH_range):
    cr_range[idx] = eos.get_c_r(PH)

critical_aspect = timing_range.critical_aspect_ratio(cr_range)

idx_crit_aspect = find_nearest(critical_aspect,1)

cm = 1/2.54
# Plot critical aspect ratio vs impact speed (Fig. 4.7)
fig,ax = plt.subplots(figsize = (30*cm,20*cm))
ax.plot(u0_range,critical_aspect)
ax.plot(u0_range[idx_crit_aspect],critical_aspect[idx_crit_aspect],'ro')
ax.plot(np.array([u0_range[idx_crit_aspect],u0_range[idx_crit_aspect]]),np.array([0,critical_aspect[idx_crit_aspect]]),'r--')
ax.fill_between(u0_range,critical_aspect,np.zeros_like(critical_aspect), color = 'red',alpha = 0.1)
ax.fill_between(u0_range,critical_aspect, plt.ylim()[1], color = 'green', alpha = 0.1)
ax.text(25,0.75,'No Steady Flow Regime', color = 'red',fontsize= 'xx-large')
ax.text(150,2.5,'Steady Flow Regime', color = 'green',fontsize= 'xx-large')
ax.grid()
ax.set_xlabel(r'Impact velocity, $u_0$ [m/s]',fontsize= 'xx-large')
ax.set_ylabel(r'Critical Aspect Ratio, $(L/D)_c$ [-]',fontsize= 'xx-large')
ax.set_xlim(xmin = 0, xmax = u0_range[-1])
ax.set_ylim(ymin = 0, ymax = critical_aspect[1])
ax.set_xlim(xmin = 0)
ax.set_ylim(ymin = 0)
plt.show()
