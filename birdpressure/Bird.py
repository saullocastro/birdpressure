import numpy as np

class Bird:
    """ The Bird class defines the material parameters as well as the dimensions.

    Notes
    -----
    All parameters are taken from Wilbeck for most accurate load prediction.

    Parameters
    ----------
    inital_density: int, default: 950
        initial density [kg/m^3]
    speed_of_sound: float, default: 1482.9
        speed of sound [m/s]
    mass: float, default: 1.814
        mass (4lb) [kg]
    aspect_ratio: float, default: 2
        length/diamter
    k: float, default: 2
        material constant for linear hugonoit

    Attributes
    ----------
    inital_density: int, default: 950
        initial density [kg/m^3]
    speed_of_sound: float, default: 1482.9
        speed of sound [m/s]
    mass: float, default: 1.814
        mass (4lb) [kg]
    aspect_ratio: float, default: 2
        length/diamter
    k: float, default: 2
        material constant for linear hugonoit
    length: float
        length of bird defined by parameters [m]
    diameter: float
        diameter of bird defined by parameters [m]

    """

    def __init__(self, initial_density=950, speed_of_sound=1482.9, mass=1.814, aspect_ratio=2, k=2):
        self.initial_density = initial_density
        self.speed_of_sound = speed_of_sound
        self.mass = mass
        self.k = k
        self.aspect_ratio = aspect_ratio
        self.length = (16 *self.mass / (np.pi*self.initial_density)) ** 1 / 3
        self.diameter = self.length / self.aspect_ratio

    def to_dict(self):
        """dict: Dictionary containing all attributes."""
        return {
            "initial_density": self.initial_density,
            "speed_of_sound": self.speed_of_sound,
            "mass": self.mass,
            "aspect_ratio": self.aspect_ratio,
            "length": self.length,
            "diameter": self.diameter
        }

    def to_string(self):
        """str: String containing all attributes."""
        return f"""
        initial_density: {self.initial_density}
        speed_of_sound: {self.speed_of_sound}
        mass: {self.mass}
        aspect_ratio: {self.aspect_ratio}
        length: {self.length}
        diameter: {self.diameter}
        """

    def get_diameter(self):
        return self.diameter

    def get_length(self):
        return self.length

    def get_mass(self):
        return self.mass

    def get_sound_speed(self):
        return self.speed_of_sound

    def get_density(self):
        return self.initial_density

    def get_LpD(self):
        return self.aspect_ratio
    def get_k(self):
        return self.k


