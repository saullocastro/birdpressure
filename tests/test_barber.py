import sys
sys.path.append('..')

import matplotlib.pyplot as plt
import numpy as np

import birdpressure.barber as barber


def test_barber():
    ######## Define Bird Impact Scenario ########
    #inputs
    theta = np.pi/2 -40*np.pi/180 # incident angle [rad]
    b = 1                         # minor axis radius [m]
    n = 39                        # grid size [-]
    U_inf = 100                   # impact speed [m/s]

    # Calculate major axis
    a = b / np.sin(theta)           # major axis radius [m]

    # Create grid/mesh
    nodes, centers, cells = barber.get_grid(a,n)

    # Define content of boundary
    cells, ids, bound_coord = barber.get_bound(centers,a,b,cells)

    # Create grid for plotting
    X, Y = np.meshgrid(centers[:,1],centers[:,0])


    ######## Plot mesh, elliptical boundary and mark cells  in boundary ########
    # Define elliptical cross-section for plot
    z_el, y_el, r_el =  barber.get_cross_ellipse(np.arange(0,2*np.pi+0.05,0.05),a,b)

    # Test index to find cell corners
    zeta_1, zeta_2, eta_1, eta_2 = barber.get_nodes(0,0,nodes)

    cm = 1/2.54
    fig,ax = plt.subplots(figsize = (50*cm,50*cm))
    ax.plot(z_el, y_el)
    ax.plot(bound_coord[:,1], bound_coord[:,0],'.')

    ax.plot(zeta_1,eta_1,'o',color= 'red')
    ax.plot(zeta_1,eta_2,'o',color= 'green')
    ax.plot(zeta_2,eta_2,'o',color= 'blue')
    ax.plot(zeta_2,eta_1,'o',color= 'black')
    ax.plot(centers[:,1][1],centers[:,0][1],'o')

    ax.set_xticks(nodes[:,1], minor=True)
    ax.set_yticks(nodes[:,0], minor=True)
    ax.set_xticks(nodes[:,1][::2])
    ax.set_yticks(nodes[:,0][::2])
    ax.grid(which = 'both', linewidth = 2)
    ax.set_xlim(a, -a)
    ax.set_ylim(ymin = -a, ymax = a)
    ax.set_xlabel('Z-Axis')
    ax.set_ylabel('Y-Axis')

    plt.savefig('test_barber_1.pdf')

    ######## Calculate pressure distribution ########
    # Determine Scource locations
    corners = barber.get_corners(ids,nodes)

    #Determine velocities
    U,V,W = barber.get_UVW (ids,centers,corners, U_inf, theta)

    #Determine Preassure Coefficents
    c_p = 1-(V**2+W**2)/(U_inf**2)

    #Assign Values to cells for Plotting
    cp_cells = cells
    cells_V = np.zeros_like(cells)
    cells_W = np.zeros_like(cells)
    for idx,id in enumerate(ids):
        cp_cells[id[0],id[1]] = c_p[idx]
        cells_V[id[0], id[1]] = V[idx]
        cells_W[id[0], id[1]] = W[idx]

    ######## Plot preassure coefficent over mesh ########
    fig, ax = plt.subplots(figsize = (30*cm,30*cm))

    ax.set_aspect('equal')
    cf = ax.pcolormesh(X,Y,cp_cells)
    fig.colorbar(cf, ax=ax)

    plt.savefig('test_barber_2.pdf')

    ######## Plot preassure coefficent over major & minor axis ########
    # Find cell ids for major & minor axis
    id_major = []
    id_minor = []

    for id in ids:
        x = 0
        y = round(centers[id[0], 0],9)
        z = round(centers[id[1], 1],9)

        if z ==0:
            id_major.append(id)
        if y ==0:
            id_minor.append(id)

    id_major = np.array(id_major)
    id_minor = np.array(id_minor)

    # Find cell cp and axis position
    cp_major = []
    axis_major = []

    cp_minor = []
    axis_minor = []

    for id in id_major:

        cp_major.append(cp_cells[id[0],id[1]])
        axis_major.append(centers[id[0],0])

    for id in id_minor:
        cp_minor.append(cp_cells[id[0], id[1]])
        axis_minor.append(centers[id[1], 1])


    cp_major = np.array(cp_major)
    axis_major = np.array(axis_major)/np.max(axis_major)

    cp_minor = np.array(cp_minor)
    axis_minor = np.array(axis_minor)/np.max(axis_minor)

    # Plotting values
    fig,ax = plt.subplots(figsize = (20*cm,20*cm))
    ax.plot(axis_major, cp_major, color = 'green', label = 'Major Axis')
    ax.plot(axis_minor, cp_minor, color = 'blue', label = 'Minor Axis')

    #ax.set_xlim(a, -a)
    #ax.set_ylim(ymin = -a, ymax = a)
    ax.set_xlabel('Nondimensional radius')
    ax.set_ylabel('C_p [-]')
    ax.legend()

    plt.savefig('test_barber_3.pdf')
