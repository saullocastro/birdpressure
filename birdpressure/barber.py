import numpy as np

# Coordinate system of impact y-z plane
#       _ _ _
#      /  ^  \
#     |   |y  |
#     | <-*   |
#      \_ z_ _/
#


#Define elliptical cross-section
##  in: angle[array/float] = angle from hor. axis, a[float] = major axis, b[float] = minor axis
## out: z[array/float] = z coordinates, y[array/float] = y coord., r[array/ float] = radius at angle
def get_cross_ellipse(angle,a,b):
    r = a*b / np.sqrt( (b * np.sin(angle))**2 + (a*np.cos(angle))**2)
    y = r * np.sin(angle)
    z = r * np.cos(angle)

    return z,y,r

# Defining Source locations/ nodes for each cell
#(eta_1,zeta_1) *_______*(eta_1,zeta_2)
#               |   ^   |
#               |   |y  |
#               | <-*   |
#(eta_2,zeta_1) *___z___* (eta_2,zeta_2)

# Defining nodes for cells/elements
##  in: i[int] = row, j[int] = column, nodes[array(n+1,2)] = coordinates of nodes (y,z)
## out: zeta_1[float] = z coord of left nodes , zeta_2[float] = z coord of right nodes ,
#       eta_1[float] = y coord of top nodes , eta_2[float] = z coord of bottom nodes
def get_nodes(i,j,nodes):
    zeta_1 = nodes[:,1][i]
    zeta_2 = nodes[:,1][i+1]

    eta_1 = nodes[:,0][j]
    eta_2 = nodes[:,0][j + 1]

    return zeta_1, zeta_2, eta_1, eta_2

# Find corners of all cells in boundary
##  in: ids[array(x,2)] = ids of cells within boundary (row, column), nodes[array(n+1,2)] = coordinates of nodes (y,z)
## out: corners[array(x,4)] = corner coordinates (eta_1,eta_2,zeta_1,zeta_2)
#       x is the number of cells/ elements within boundary
def get_corners(ids,nodes):

    corners = np.zeros([ids.shape[0],4])

    for n,id in enumerate(ids):

        zeta_1 = nodes[:, 1][id[1]]
        zeta_2 = nodes[:, 1][id[1] + 1]

        eta_1 = nodes[:, 0][id[0]]
        eta_2 = nodes[:, 0][id[0] + 1]
        corners[n,:] = [eta_1, eta_2, zeta_1, zeta_2]
    return corners

# Calculating velocity field
##  in: x,y,z[float] = coordinates of point of interest, eta_1,eta_2,zeta_1,zeta_2[float] = source coordinates
## out: u,v,w[float] = velocity contribution in x,y,z,
#       r_1,r_2,r_3,r_4[float] = distance from point of interest to source locations
def get_velocity(x,y,z,eta_1, eta_2, zeta_1, zeta_2):

    r_1 = np.sqrt(x**2 + (y - eta_1)**2 + (z - zeta_1)**2)
    r_2 = np.sqrt(x**2 + (y - eta_2)**2 + (z - zeta_1)**2)
    r_3 = np.sqrt(x**2 + (y - eta_2)**2 + (z - zeta_2)**2)
    r_4 = np.sqrt(x**2 + (y - eta_1)**2 + (z - zeta_2)**2)

    if x ==0 :
        u = (np.pi/2*(((z - zeta_2) * (y - eta_2))/abs((z - zeta_2) * (y - eta_2)))  + \
        np.pi/2*((z - zeta_1) * (y - eta_1)) / abs((z - zeta_1) * (y - eta_1)) - \
        np.pi/2*((z - zeta_1) * (y - eta_2) / abs((z - zeta_1) * (y - eta_2))) - \
        np.pi/2*((z - zeta_2) * (y - eta_1) / abs((z - zeta_2) * (y - eta_1))))
    else:

        u = (np.arctan((z-zeta_2) * (y-eta_2) / (x*r_3)) + np.arctan((z-zeta_1) * (y-eta_1) / (x*r_1)) -\
        np.arctan((z-zeta_1) * (y-eta_2) / (x*r_2)) - np.arctan((z-zeta_2) * (y-eta_1) / (x*r_4)))

    v =  np.log( ((r_3+(zeta_2-z))*(r_1+(zeta_1-z))) / ((r_4+(zeta_2-z))*(r_2+(zeta_1-z))))
    w =  np.log( ((r_3+( eta_2-y))*(r_1+( eta_1-y))) / ((r_2+( eta_2-y))*(r_4+( eta_1-y))))


    return u, v, w, r_1,r_2,r_3,r_4

# Define grid/mesh
##  in: a[float] = major axis, n[int] = number of elemts per row/column
## out: nodes[array(n+1,2)] = coordinates of nodes (y,z), centers[array(n,2)] = coordinates of centers (z,y),
#       cells[array(n,n)] = array representing each cell/element in grid
def get_grid(a,n):
    z_nodes = np.linspace(-a, a, n + 1)
    y_nodes = np.linspace(a, -a, n + 1)
    nodes = np.transpose(np.vstack([y_nodes, z_nodes]))

    cells = np.zeros([n, n])
    cell_size = abs(z_nodes[1] - z_nodes[0])

    center_z = np.linspace(-a + cell_size / 2, a - cell_size / 2, n)
    center_y = np.linspace(a - cell_size / 2, -a + cell_size / 2, n)
    centers = np.transpose(np.vstack([center_y, center_z]))

    return nodes, centers, cells

#Defining which cells/elemts are within elliptical boundary
##  in: centers[array(n,2)] = coordinates of centers (z,y), a[float] = major axis, b[float] = minor axis,
#       cells[array(n,n)] = array representing each cell/element in grid
## out: cells[array(n,n)] = cells in boundary are value 1, ids[array(x,2)] = ids of cells within boundary (row, column),
#       bound_coord[array(x,2)] = coordinates of cell centers within boundary (z,y),
#       x is the number of cells/ elements within boundary
def get_bound(centers,a,b,cells):

    n = cells.shape[0]

    z_coord = []
    y_coord = []
    ids = []

    for idx_y in np.arange(n):
        for idx_z in np.arange(n):


            if (n-1)%2==0 and idx_z == (n-1)/2 and idx_y ==(n-1)/2:
                print('hallo')
                cells[idx_z, idx_y] = 1
                z_coord.append(centers[:, 1][idx_z])
                y_coord.append(centers[:, 0][idx_y])
                ids.append([idx_y, idx_z])

            else:

                r = np.sqrt(centers[:, 1][idx_z] ** 2 + centers[:, 0][idx_y] ** 2)
                if centers[:, 1][idx_z] ==0:
                    angle = np.pi/2 * centers[:, 0][idx_y]/abs(centers[:, 0][idx_y])
                else:
                    angle = np.arctan(centers[:, 0][idx_y] / centers[:, 1][idx_z])
                z_bound, y_bound, r_bound = get_cross_ellipse(angle, a, b)

                if r <= r_bound:
                    cells[idx_y, idx_z] = 1
                    z_coord.append(centers[:, 1][idx_z])
                    y_coord.append(centers[:, 0][idx_y])
                    ids.append([idx_y, idx_z])



    ids = np.array(ids)
    bound_coord = np.transpose(np.vstack([y_coord, z_coord]))

    return cells, ids, bound_coord

# Determine Velocities over y-z plane
##  in: ids[array(x,2)] = ids of cells within boundary (row, column), centers[array(n,2)] = coordinates of centers (z,y),
#       corners[array(x,4)] = corner coordinates (eta_1,eta_2,zeta_1,zeta_2), U_inf[float] = Impact Velocity
#       theta[float] = Impact angle
## out: U[array(x)] = array of velocities in x direction (out of y-z plane), V[array(x)] = array of velocities in y direction,
#       W[array(x)] = array of velocities in z direction
#       x is the number of cells/ elements within boundary
def get_UVW (ids,centers,corners, U_inf, theta):
    uvw = np.zeros([ids.shape[0], 3])
    for idx_1, id_cent in enumerate(ids):

        x = 0
        y = centers[id_cent[0], 0]
        z = centers[id_cent[1], 1]

        for idx in np.arange(corners.shape[0]):
            eta_1 = corners[idx, 0]
            eta_2 = corners[idx, 1]
            zeta_1 = corners[idx, 2]
            zeta_2 = corners[idx, 3]

            u, v, w, r_1, r_2, r_3, r_4 = get_velocity(x, y, z, eta_1, eta_2, zeta_1, zeta_2)
            uvw[idx_1, :] += [u, v, w]

    U = U_inf * np.sin(theta) + U_inf * np.sin(theta) / (2 * np.pi) * uvw[:, 0]
    V = U_inf * np.cos(theta) + U_inf * np.sin(theta) / (2 * np.pi) * uvw[:, 1]
    W = U_inf * np.sin(theta) / (2 * np.pi) * uvw[:, 2]

    return U,V,W