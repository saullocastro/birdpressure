from birdpressure.Bird import Bird
from birdpressure.ImpactScenario import ImpactScenario
from birdpressure.Timing import Timing

import numpy as np
import matplotlib.pyplot as plt

u0_range = np.linspace(0.1,1500.1,41)
bird = Bird()
scen = ImpactScenario(u0_range, 90, use_radian = False)
tim = Timing(bird,scen)

P_us = tim.det_peak_P()
P_u0 = tim.det_peak_P(use_us = False)

cm = 1/2.54

fig,ax = plt.subplots(figsize = (30*cm,20*cm))
ax.plot(u0_range,P_u0 * 10**-6, label = r'$\rho*c_0*u_0$')
ax.plot(u0_range,P_us * 10**-6, label = r'$\rho*u_s*u_0$')
ax.legend(fontsize= 'xx-large', loc = 'upper left')
ax.grid()
ax.set_xlabel(r'Impact velocity, $u_0$ [m/s]',fontsize= 'xx-large')
ax.set_ylabel(r'Hugoniot pressur, $P_H$ [MPa]',fontsize= 'xx-large')
ax.set_xlim(xmin = 0)
ax.set_ylim(ymin = 0)
plt.show()