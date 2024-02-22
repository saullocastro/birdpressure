import numpy as np

class Bird:

    # make docu
    " Bird class for the generation of the bird object"

    def __init__(self, initial_density=950, speed_of_sound=1482.9, mass=1.814, aspect_ratio=2, k=2):
        self.initial_density = initial_density
        self.speed_of_sound = speed_of_sound
        self.mass = mass
        self.k = k
        self.aspect_ratio = aspect_ratio
        self.length = (16 / np.pi * self.mass / self.initial_density) ** 1 / 3
        self.diameter = self.length / self.aspect_ratio

    def to_dict(self):
        return {
            "initial_density": self.initial_density,
            "speed_of_sound": self.speed_of_sound,
            "mass": self.mass,
            "aspect_ratio": self.aspect_ratio,
            "length": self.length,
            "diameter": self.diameter
        }

    def to_string(self):
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
        return self.m

    def get_sound_speed(self):
        return self.speed_of_sound

    def get_density(self):
        return self.initial_density

    def get_LpD(self):
        return self.aspect_ratio
    def get_k(self):
        return self.k


