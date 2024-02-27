from birdpressure.Bird import Bird
from birdpressure.ImpactScenario import ImpactScenario
from birdpressure.Timing import Timing

import matplotlib.pyplot as plt
import numpy as np

def test_oblique_impact():
    u0_range = np.linspace(0.1,1500.1,101)
    angle_range = np.linspace(15,90,4)
    bird = Bird()


    scen_15 = ImpactScenario(u0_range, 15, use_radian = False)
    scen_30 = ImpactScenario(u0_range, 30, use_radian = False)
    scen_45 = ImpactScenario(u0_range, 45, use_radian = False)
    scen_90 = ImpactScenario(u0_range, 90, use_radian = False)

    tim_15 = Timing(bird,scen_15)
    PH_15  = tim_15.det_peak_P()

    tim_30 = Timing(bird,scen_30)
    PH_30  = tim_30.det_peak_P()

    tim_45 = Timing(bird,scen_45)
    PH_45  = tim_45.det_peak_P()

    tim_90 = Timing(bird,scen_90)
    PH_90  = tim_90.det_peak_P()

    cm = 1/2.54
    fig,ax = plt.subplots(figsize = (30*cm,20*cm))

    ax.plot(u0_range, PH_15*10**(-6), label= r'$\alpha$ = $15^{\circ}$')
    ax.plot(u0_range, PH_30*10**(-6), label= r'$\alpha$ = $30^{\circ}$')
    ax.plot(u0_range, PH_45*10**(-6), label= r'$\alpha$ = $45^{\circ}$')
    ax.plot(u0_range, PH_90*10**(-6), label= r'$\alpha$ = $90^{\circ}$')
    ax.grid()
    ax.legend(fontsize= 'xx-large', loc = 'upper left')
    ax.set_xlabel(r'Impact Velocity, $u_0$ [m/s]',fontsize= 'xx-large')
    ax.set_ylabel(r'Hugonoit Pressure, $P_H$ [MPa]',fontsize= 'xx-large')
    ax.set_xlim(xmin = 0, xmax = u0_range[-1])
    ax.set_ylim(ymin = 0)
    plt.show()