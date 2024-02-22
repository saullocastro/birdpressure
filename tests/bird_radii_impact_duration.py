from BirdPressure.Bird import Bird
from BirdPressure.ImpactScenario import ImpactScenario
from BirdPressure.Timing import Timing
from BirdPressure.EOS import EOS

import numpy as np
import matplotlib.pyplot as plt


u0_range = np.linspace(0.1,300.1,41)
bird = Bird()
bird.diameter = 2 * np.array([0.01,0.02,0.03])
scen = ImpactScenario(u0_range, 90, use_radian = False)
timing_range = Timing(bird,scen)
PH_range = timing_range.det_peak_P()
EOS_pol = EOS(bird, 'Polynomial_EOS')

pol_cr_range = np.zeros_like(PH_range)
tb_1 = np.zeros_like(PH_range)
tb_2 = np.zeros_like(PH_range)
tb_3 = np.zeros_like(PH_range)

for idx,PH in enumerate(PH_range):
    pol_cr_range[idx] = EOS_pol.get_c_r(PH)
    tb_range = timing_range.find_t_b(pol_cr_range[idx])
    tb_1[idx] = tb_range[0]
    tb_2[idx] = tb_range[1]
    tb_3[idx] = tb_range[2]


cm = 1/2.54
fig,ax = plt.subplots(figsize = (30*cm,20*cm))
ax.plot(u0_range,  tb_1*10e5, label='a = 0.01')
ax.plot(u0_range,  tb_2*10e5, label='a = 0.02')
ax.plot(u0_range,  tb_3*10e5, label='a = 0.03')
ax.legend(fontsize= 'xx-large', loc = 'upper right')
ax.grid()
ax.set_xlabel(r'Impact velocity, $u_0$ [m/s]',fontsize= 'xx-large')
ax.set_ylabel(r'Impact Duaration, [$\mu$ s]',fontsize= 'xx-large')
ax.set_xlim(xmin = 0)
ax.set_ylim(ymin = 0)
plt.show()




