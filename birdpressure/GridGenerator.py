from birdpressure.Bird import Bird
from birdpressure.ImpactScenario import ImpactScenario
import matplotlib.pyplot as plt
import numpy as np

class GridGenerator:
    r"""Defines grid and boundary of initial impact area

    Note
    ----
    Coordinate system of impact y-z plane:
                      _ _ _
                     /  ^  \
                    |   |y  |
                    | <-*   |
                    \_ z_ _/
    Defining Elements/Cells:
    (eta_1,zeta_1) *_______*(eta_1,zeta_2)
                   |   ^   |
                   |   |y  |
                   | <-*   |
    (eta_2,zeta_1) *___z___* (eta_2,zeta_2)

    Parameters
    ----------
    bird: class
    impact_scenario: class
    n_elements: int
        number of elements per row/column

    Attributes
    ----------
    n_elements: int
        number of elements per row/column
    maj_axis: float
        half the length of elliptical impact area major axis [m]
    min_axis: float
        half the length of elliptical impact area major axis [m]
    nodes: ndarray
         array containing all coordinates of nodes [y,z] (n_elements+1,2)
    centers: ndarray
        array containing all coordinates of element centers [y,z] (n_elements,2)
    cells_bound: ndarray
        array representing grid with 1 values for cell/element in bounds (n_elements, n_elements)
    ids: ndarray
        array containing element id of all elements withing impact boundary [row, column]
    corners: ndarray
        array containing all corner coordinates of elements in boundary [eta_1,eta_2,zeta_1,zeta_2]
    resolution: float
        elements size
    mesh_X: ndarray
        array created by np.mesh() with x coordinates (n_elements, n_elements)
    mesh_Y: ndarray
        array created by np.mesh() with y coordinates (n_elements, n_elements)
    impact_area: float
        numerically determined impact area [m^2]
    exact_impact_area : float
        exact impact area [m^2]

    Other Parameters
    ----------------
    cells: ndarray
        array of zeros representing each cell/element in grid (n_elements, n_elements)
    bound_coord: ndarray
        coordinates of cell centers within boundary [z,y]
    """
    def __init__(self, bird: Bird, impact_scenario: ImpactScenario, n_elements: int):

        self.n_elements = n_elements
        self.min_axis = bird.get_diameter() / 2
        self.maj_axis = self.min_axis / np.sin(impact_scenario.get_impact_angle())
        self.nodes, cells, self.centers = self.get_grid()
        self.cells_bound, self.ids, bound_coord = self.get_bound(cells)
        self.corners = self.get_corners(self.nodes)
        self.resolution =  2*self.maj_axis/self.n_elements
        self.mesh_X, self.mesh_Y = np.meshgrid(self.centers[:, 1], self.centers[:, 0])
        self.impact_area = self.ids.shape[0] * self.resolution**2
        self.exact_impact_area = np.pi * self.min_axis * self.maj_axis
    def get_grid(self):
        """ This function generates the Grid

        Returns
        -------
        nodes: ndarray
            array containing all coordinates of nodes [y,z] (n_elements+1,2)
        cells: ndarray
            array representing each cell/element in grid (n_elements, n_elements)
        centers: ndarray
            array containing all coordinates of element centers [y,z] (n_elements,2)
        """

        z_nodes = np.linspace(-self.maj_axis, self.maj_axis, self.n_elements + 1)
        y_nodes = np.linspace(self.maj_axis, -self.maj_axis, self.n_elements + 1)
        nodes = np.transpose(np.vstack([y_nodes, z_nodes]))

        cells = np.zeros([self.n_elements, self.n_elements])
        cell_size = abs(z_nodes[1] - z_nodes[0])

        center_z = np.linspace(-self.maj_axis + cell_size / 2, self.maj_axis - cell_size / 2, self.n_elements)
        center_y = np.linspace(self.maj_axis - cell_size / 2, -self.maj_axis + cell_size / 2, self.n_elements)
        centers = np.transpose(np.vstack([center_y, center_z]))
        return nodes, cells, centers


    def get_bound(self, cells):
        """ This function determines which elements/cells are within impact area

        Parameters
        ----------
        cells: ndarray
            array representing each cell/element in grid (n_elements, n_elements)

        Returns
        -------
        cells_bound: ndarray
            array representing grid with 1 values for cell/element in bounds (n_elements, n_elements)
        ids: ndarray
            array containing element id of all elements withing impact boundary [row, column]
        bound_coord: ndarray
            coordinates of cell centers within boundary [z,y]
        """

        z_coord = []
        y_coord = []
        ids = []

        for idx_y in np.arange(self.n_elements):
            for idx_z in np.arange(self.n_elements):

                if (self.n_elements - 1) % 2 == 0 and idx_z == (self.n_elements - 1) / 2 and idx_y == (
                        self.n_elements - 1) / 2:
                    cells[idx_z, idx_y] = 1
                    z_coord.append(self.centers[:, 1][idx_z])
                    y_coord.append(self.centers[:, 0][idx_y])
                    ids.append([idx_y, idx_z])

                else:

                    r = np.sqrt(self.centers[:, 1][idx_z] ** 2 + self.centers[:, 0][idx_y] ** 2)
                    if self.centers[:, 1][idx_z] == 0:
                        angle = np.pi / 2 * self.centers[:, 0][idx_y] / abs(self.centers[:, 0][idx_y])
                    else:
                        angle = np.arctan(self.centers[:, 0][idx_y] / self.centers[:, 1][idx_z])
                    z_bound, y_bound, r_bound = self.get_cross_ellipse(angle)

                    if r <= r_bound:
                        cells[idx_y, idx_z] = 1
                        z_coord.append(self.centers[:, 1][idx_z])
                        y_coord.append(self.centers[:, 0][idx_y])
                        ids.append([idx_y, idx_z])

        ids = np.array(ids)
        bound_coord = np.transpose(np.vstack([y_coord, z_coord]))
        cells_bound = cells

        return cells_bound, ids, bound_coord


    def get_cross_ellipse(self, angle):
        """This function creates coordinates of impact area boundary for given angle from minor axis

        Parameters
        ----------
        angle: float or ndarray
            angle from minor axis [rad]

        Returns
        -------
        z: float or ndarray
            z coord
        y: float or ndarray
            y coord
        r: float or ndarray
            distance from center
        """
        r = self.maj_axis * self.min_axis / np.sqrt(
            (self.min_axis * np.sin(angle)) ** 2 + (self.maj_axis * np.cos(angle)) ** 2)
        y = r * np.sin(angle)
        z = r * np.cos(angle)

        return z, y, r


    def get_corners(self, nodes: np.ndarray):
        """Determines corner coord for given element/cell id

        Parameters
        ----------
        nodes: nd.array
            element/cell id [row, column]

        Returns
        -------
        corners: ndarray
            array containing all corner coordinates of elements in boundary [eta_1,eta_2,zeta_1,zeta_2]
        """
        corners = np.zeros([self.ids.shape[0], 4])

        for n, id in enumerate(self.ids):
            zeta_1 = nodes[:, 1][id[1]]
            zeta_2 = nodes[:, 1][id[1] + 1]

            eta_1 = nodes[:, 0][id[0]]
            eta_2 = nodes[:, 0][id[0] + 1]
            corners[n, :] = [eta_1, eta_2, zeta_1, zeta_2]

        return corners

    def show_grid(self):
        """Function that plots grid and impact area."""
        z_el, y_el, r_el = self.get_cross_ellipse(np.linspace(0,2*np.pi,30))


        cm = 1 / 2.54
        fig, ax = plt.subplots(figsize=(20 * cm, 20 * cm))
        ax.plot(z_el, y_el,linewidth = 3, label = 'Boundary')

        ax.set_xticks(self.nodes[:, 1], minor=True)
        ax.set_yticks(self.nodes[:, 0], minor=True)
        ax.set_xticks(self.nodes[:, 1][::2])
        ax.set_yticks(self.nodes[:, 0][::2])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.grid(which='both', linewidth=2)
        ax.set_xlim(self.maj_axis, -self.maj_axis)
        ax.set_ylim(ymin=-self.maj_axis, ymax=self.maj_axis)
        ax.set_xlabel('Z-Axis')
        ax.set_ylabel('Y-Axis')
        ax.legend()

        plt.show()

