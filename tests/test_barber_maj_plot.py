import sys
sys.path.append('..')

import matplotlib.pyplot as plt
import numpy as np

import birdpressure.barber as barber


def test_barber_maj_plot():
    ######## Define Bird Impact Scenario ########
    #inputs
    thetas = np.array([25,45,60,75,85]) * np.pi/180 # incident angles [rad]
    b = 1                                           # minor axis radius [m]
    n = 79                                          # grid size [-]
    U_inf = 100                                     # impact speed [m/s]

    #Alocate memory for C_p data in form of a np.array
    cp_s_maj = np.zeros([n,len(thetas)])

    ######## Determine C_p distribution for all impact angles theta ########
    #Loop through impact angles
    for idx,theta in enumerate(thetas):
        #Define mesh & boundaries
        a                       = b / np.sin(theta)           # major axis radius [m]
        nodes, centers, cells   = barber.get_grid(a,n)
        cells, ids, bound_coord = barber.get_bound(centers,a,b,cells)
        X, Y                    = np.meshgrid(centers[:,1],centers[:,0])
        id_major                = []
    # Find ids for cells//elements along major axis
        for id in ids:
            x = 0
            y = round(centers[id[0], 0],9)
            z = round(centers[id[1], 1],9)
            if z ==0:
                id_major.append(id)
    # Use barber to calculate C_p dits.
        id_major        = np.array(id_major)
        corners         = barber.get_corners(ids,nodes)
        U,V,W           = barber.get_UVW (id_major,centers,corners, U_inf, theta)
        cp_maj          = 1-(V**2+W**2)/(U_inf**2)
        cp_s_maj[:,idx] = cp_maj


    ######## Plot C_p distribution for changing impact angles ########
    #Normalise major axis
    axis_maj = np.linspace(-1,1,len(cp_maj))
    #Convert inch to cm
    cm      = 1/2.54

    #Plot
    fig,ax  = plt.subplots(figsize = (50*cm,50*cm))

    for idx,theta in enumerate(thetas):
        ax.plot(axis_maj, cp_s_maj[:,idx], label = str(round(theta * 180/np.pi)) + '[deg]')

    ax.set_title('Pressure coefficient  vs. nondimensional radius along the major axis')
    ax.set_xlabel('Nondimensional radius')
    ax.set_ylabel('C_p [-]')
    ax.legend()

    plt.savefig('test_barber_maj_plot_1.pdf')
