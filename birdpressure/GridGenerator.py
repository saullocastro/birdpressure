from BirdPressure.Bird import Bird
from BirdPressure.ImpactScenario import ImpactScenario
import matplotlib.pyplot as plt
import numpy as np

class GridGenerator:
    def __init__(self, bird: Bird, impact_scenario: ImpactScenario, n_elements: int):
        """
        Grid Generator for Bird Impact Simulation
        :param bird: Bird Class
        :param impact_scenario: ImpactScenario Class
        :param n_elements: number of elemts per row/column
        """
        self.n_elements = n_elements
        self.min_axis = bird.get_diameter() / 2
        self.maj_axis = self.min_axis / np.sin(impact_scenario.get_impact_angle())
        self.nodes, cells, self.centers = self.get_grid()
        cells, self.ids, bound_coord = self.get_bound(cells)
        self.corners = self.get_corners(self.nodes)
        self.resolution =  2*self.maj_axis/self.n_elements
        self.mesh_X, self.mesh_Y = np.meshgrid(self.centers[:, 1], self.centers[:, 0])
        self.impact_area = self.ids.shape[0] * self.resolution**2
    def get_grid(self):
        """
                This function generates the Grid

                Params:
                    - maj_axis (float): this is an angle
                    - n_elements (integer) : number of elemts per row/column

                Return:
                    - node (numpy array [n_elements+1, 2]) : coordinates of nodes [y,z]
                    - cells (numpy array [n_elements, n_elements]): array representing each cell/element in grid
                    - centers (numpy array [n_elements, 2]): coordinates of centers [y,z]

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
        """

        Parameters
        ----------
        cells

        Returns
        -------

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

        return cells, ids, bound_coord


    def get_cross_ellipse(self, angle):
        """"
        Hallo I am Jonas

        Params:
            - angle (float): this is an angle

        Return:
            - lala

        """"
        r = self.maj_axis * self.min_axis / np.sqrt(
            (self.min_axis * np.sin(angle)) ** 2 + (self.maj_axis * np.cos(angle)) ** 2)
        y = r * np.sin(angle)
        z = r * np.cos(angle)

        return z, y, r


    def get_corners(self, nodes: np.ndarray):
        """

        Parameters
        ----------
        nodes

        Returns
        -------

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
        """

        """
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

