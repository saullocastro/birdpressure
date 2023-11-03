import sys
sys.path.append('..')

import numpy as np
import matplotlib.pyplot as plt

import birdpressure.wilbeck as wilbeck


def test_wilbeck():
    ######## Define Bird ########
    # Material Prop. taken from Wilbeck:

    k   = 2       #[-]
    c_0 = 1482.9  #[m/s]
    rho = 950     #[kg/m^3]

    # Bird Dimension
    L   = 2                  #[m]
    D   = 1                  #[m]
    M   = np.pi/4*D**2*L*rho #[kg]


    ######## Showcase effect of compressibility on P_H (Hedayati Fig. 4.3) ########
    # Set velocity range
    u_0 = np.arange(0,300,20) #[m/s]

    # Calculated uncompressible P_H
    P_H1 = rho * c_0 * u_0 #[Pa]

    # Calculated shock velocity of compressible material
    u_s = wilbeck.find_u_s(c_0,k,u_0) #[m/s]

    # Calculated compressible P_H
    P_H2 = rho * u_s * u_0 #[Pa]

    # Set inch to cm conversion
    cm = 1/2.54

    # Plot effect on P_H (Hedayati Fig. 4.3)
    fig,ax = plt.subplots(figsize = (30*cm,20*cm))
    ax.plot(u_0,P_H1 * 10**-6, label = r'$\rho*c_0*u_0$')
    ax.plot(u_0,P_H2 * 10**-6, label = r'$\rho*u_s*u_0$')
    ax.legend(fontsize= 'xx-large', loc = 'upper left')
    ax.grid()
    ax.set_xlabel(r'Impact velocity, $u_0$ [m/s]',fontsize= 'xx-large')
    ax.set_ylabel(r'Hugoniot pressur, $P_H$ [MPa]',fontsize= 'xx-large')
    ax.set_xlim(xmin = 0)
    ax.set_ylim(ymin = 0)

    plt.savefig('test_wilbeck_1.pdf')


    ######## Shock regime (Hedayati Sec. 4.2.1) ########
    # Calculated speed of sound in shocked region
    c_r = wilbeck.find_c_r(c_0,k,u_0) #[m/s]

    # Plot shock velocity vs. speed of sound (Fig. 4.5)
    fig,ax = plt.subplots(figsize = (30*cm,20*cm))
    ax.plot(u_0[0:16],u_s[0:16], label = r'$u_s$')
    ax.plot(u_0[0:16],c_r[0:16], label = r'$c_r$')
    ax.legend(fontsize= 'xx-large', loc = 'upper left')
    ax.grid()
    ax.set_xlabel(r'Impact velocity, $u_0$ [m/s]',fontsize= 'xx-large')
    ax.set_ylabel('Wave speed,[m/s]',fontsize= 'xx-large')
    ax.set_xlim(xmin = 0)
    ax.set_ylim(ymin = 0)

    plt.savefig('test_wilbeck_2.pdf')


    # Set diameter range
    a = np.array([0.01, 0.02, 0.03]) #[m]

    # Set velocity range
    u_0 = np.arange(0,360,10) #[m/s]

    # Calculated speed of sound in shocked region
    c_r= wilbeck.find_c_r(c_0,k, u_0) #[m//s]

    # PLot affect of bird diameter on impact duration (Fig. 4.6)
    fig,ax = plt.subplots(figsize = (30*cm,20*cm))
    for a_val in a:
        t_b = wilbeck.find_t_b(c_r, a_val) #[sec]
        ax.plot(u_0, t_b*10e5, label='a = ' + str(a_val))
    ax.legend(fontsize= 'xx-large', loc = 'upper right')
    ax.grid()
    ax.set_xlabel(r'Impact velocity, $u_0$ [m/s]',fontsize= 'xx-large')
    ax.set_ylabel(r'Impact Duaration, [$\mu$ s]',fontsize= 'xx-large')
    ax.set_xlim(xmin = 0)
    ax.set_ylim(ymin = 0)

    plt.savefig('test_wilbeck_3.pdf')


    ######## Release regime (Hedayati Sec. 4.2.2) ########
    # Calculated speed of sound in shocked region
    cr = wilbeck.find_c_r(c_0,k ,u_0) #[m/s]
    u_s = wilbeck.find_u_s(c_0,k,u_0) #[m/s]

    # Calculated critical bird aspect ratio
    L_crit = wilbeck.find_LD_crit(u_s,cr,u_0)

    # Find point where (L/D)_c = 1
    idx_LD = wilbeck.find_nearest(L_crit,1)

    # Plot critical aspect ratio vs impact speed (Fig. 4.7)
    fig,ax = plt.subplots(figsize = (30*cm,20*cm))
    ax.plot(u_0,L_crit)
    ax.plot(u_0[idx_LD],L_crit[idx_LD],'ro')
    ax.plot(np.array([u_0[idx_LD],u_0[idx_LD]]),np.array([0,L_crit[idx_LD]]),'r--')
    ax.fill_between(u_0,L_crit,np.zeros_like(L_crit), color = 'red',alpha = 0.1)
    ax.fill_between(u_0,L_crit, plt.ylim()[1], color = 'green', alpha = 0.1)
    ax.text(25,0.75,'No Steady Flow Regime', color = 'red',fontsize= 'xx-large')
    ax.text(150,2.5,'Steady Flow Regime', color = 'green',fontsize= 'xx-large')
    ax.grid()
    ax.set_xlabel(r'Impact velocity, $u_0$ [m/s]',fontsize= 'xx-large')
    ax.set_ylabel(r'Critical Aspect Ratio, $(L/D)_c$ [-]',fontsize= 'xx-large')
    ax.set_xlim(xmin = 0, xmax = u_0[-1])
    ax.set_ylim(ymin = 0, ymax = L_crit[1])
    ax.set_xlim(xmin = 0)
    ax.set_ylim(ymin = 0)

    plt.savefig('test_wilbeck_4.pdf')


    ######## Steady flow regime (Hedayati Sec. 4.2.3) ########
    # Set impact scenario
    a_max = 1   #[m]
    u_0 = 100   #[m/s]
    rho_0 = 950 #[m/s]

    # Create array for radius from center
    r = np.linspace(0,a_max,101) #[m]

    # Create arrays for pressure and density along radius Bank formulation
    rho_B = np.zeros_like(r)
    P_B = np.zeros_like(r)

    # Create arrays for pressure and density along radius Leach formulation
    rho_L = np.zeros_like(r)
    P_L = np.zeros_like(r)

    # Fill arrays with values
    for idx,r_val in enumerate(r):
        if idx == 0:
            P_start     = rho_0 * c_0 * u_0
            rho_start   = rho_0

        else:
            P_start     = P_B[idx-1]
            rho_start   = rho_B[idx-1]

        # Banks and Leach Formulation Solver
        result_B    = wilbeck.solve_Banks(u_0, r_val,P_start, rho_start)
        result_L    = wilbeck.solve_Leach(u_0, r_val, P_start, rho_start)

        rho_B[idx]  = result_B[0]
        P_B[idx]    = result_B[1]

        rho_L[idx]  = result_L[0]
        P_L[idx]    = result_L[1]

    #Plot Pressure distributions Banks and Chandrasekhara & Leach and Walker
    fig,ax = plt.subplots(figsize = (30*cm,20*cm))
    ax.plot(r,P_B*10**(-6),label='Banks and Chandrasekhara')
    ax.plot(r,P_L*10**(-6),label='Leach and Walker')
    ax.grid()
    ax.legend(fontsize= 'xx-large', loc = 'lower left')
    ax.set_xlabel(r'Distance from Impact center, $r$ [m]',fontsize= 'xx-large')
    ax.set_ylabel(r'Pressure, $P$ [MPa]',fontsize= 'xx-large')
    ax.set_xlim(xmin = 0,xmax = a_max)
    ax.set_ylim(ymin = 0)

    plt.savefig('test_wilbeck_5.pdf')

    #Plot pressure ditribution and density
    fig, ax1 = plt.subplots(figsize = (30*cm,20*cm))

    color = 'tab:red'
    ax1.set_xlabel(r'Distance from Impact center, $r$ [m]',fontsize= 'xx-large')
    ax1.grid()
    ax1.set_ylabel(r'Pressure, $P$ [MPa]',fontsize= 'xx-large', color = color)
    ax1.plot(r,P_B*10**(-6), color= 'indianred', label='Banks and Chandrasekhara')
    ax1.plot(r,P_L*10**(-6), color= 'maroon', label='Leach and Walker')
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(r'Density, $\rho$ [kg/$m^3$]',fontsize= 'xx-large', color=color)  # we already handled the x-label with ax1
    ax2.plot(r, rho_B, color='royalblue', label='Banks and Chandrasekhara')
    ax2.plot(r, rho_L, color='navy', label='Leach and Walker')
    ax2.legend(fontsize= 'xx-large', loc = 'upper right')
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.savefig('test_wilbeck_6.pdf')

    # Compute steady state pressure with constant rho
    P_B_const = wilbeck.find_banks(rho, u_0, 0.5, r, 0.5)
    P_L_const = wilbeck.find_leach(rho, u_0, 2.58, r, 0.5)

    #Plot pressure ditribution and density with const. rho
    fig, ax1 = plt.subplots(figsize = (30*cm,20*cm))

    color = 'tab:red'
    ax1.set_xlabel(r'Distance from Impact center, $r$ [m]',fontsize= 'xx-large')
    ax1.grid()
    ax1.set_ylabel(r'Pressure, $P$ [MPa]',fontsize= 'xx-large', color = color)
    ax1.plot(r,P_B_const*10**(-6), color= 'indianred', label='Banks and Chandrasekhara')
    ax1.plot(r,P_L_const*10**(-6), color= 'maroon', label='Leach and Walker')
    ax1.legend(fontsize= 'xx-large', loc = 'center right')
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel(r'Density, $\rho$ [kg/$m^3$]',fontsize= 'xx-large', color=color)  # we already handled the x-label with ax1
    ax2.plot(r, rho_0 * np.ones_like(P_B_const), color='royalblue', label='Banks and Chandrasekhara')
    ax2.plot(r, rho_0* np.ones_like(P_L_const), color='navy', label='Leach and Walker')
    ax2.legend(fontsize= 'xx-large', loc = 'upper right')
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.savefig('test_wilbeck_7.pdf')

    ######## Oblique impact (Hedayati Sec. 4.3.2) ########
    # Set density
    rho = 950     #[kg/m^3]

    # Set impact scenarios
    alpha_val   = np.array([90,45,30,15]) #[deg]
    u_0         = np.arange(0,300,20)     #[m/s]

    # Plot P_H for different impact angles (Fig. 4.1)
    fig,ax = plt.subplots(figsize = (30*cm,20*cm))

    for alpha in alpha_val:
        u_ob = u_0 * np.sin(alpha * np.pi/180)
        u_s = wilbeck.find_u_s(c_0,k,u_ob)
        P_H = rho * u_s * u_ob
        ax.plot(u_0, P_H*10**(-6), label= r'$\alpha$ =' + str(alpha) + r'$^{\circ}$')

    ax.grid()
    ax.legend(fontsize= 'xx-large', loc = 'upper left')
    ax.set_xlabel(r'Impact Velocity, $u_0$ [m/s]',fontsize= 'xx-large')
    ax.set_ylabel(r'Hugonoit Pressure, $P_H$ [MPa]',fontsize= 'xx-large')
    ax.set_xlim(xmin = 0, xmax = u_0[-1])
    ax.set_ylim(ymin = 0)

    plt.savefig('test_wilbeck_8.pdf')

