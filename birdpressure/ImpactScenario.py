import numpy as np


class ImpactScenario:
    """This class defines the impact scenario

    Parameters
    ----------
    impact_velocity: float or ndarray
        impact velocity [m/s]
    impact_angle: float or ndarray
        impact angle [rad/deg]
    use_radian: bool, default: True
        indicates whether you use radian or degree

    Attributes
    ----------
    impact_velocity: float or ndarray
        impact velocity [m/s]
    impact_angle: float or ndarray
        impact angle [rad/deg]
    use_radian: bool, default: True
        indicates whether you use radian or degree


    """
    def __init__(self, impact_velocity, impact_angle, use_radian = True):
        self.impact_velocity = impact_velocity
        self.use_radian = use_radian

        if use_radian:
            self.impact_angle = impact_angle
        else:
            self.impact_angle = impact_angle * np.pi/180

        self.normal_velocity = self.impact_velocity * np.sin(self.impact_angle)

    def get_impact_velocity(self):
        return self.impact_velocity

    def get_normal_velocity(self):
        return self.normal_velocity

    def get_impact_angle(self):
        return self.impact_angle

    def to_dict(self):
        """dict: Dictionary containing all attributes."""

        return {
            "impact_velocity": self.impact_velocity,
            "impact_angle": self.impact_angle,
            "normal_velocity": self.normal_velocity
        }

    def to_string(self):
        """str: String containing all attributes."""
        return f"""
        impact_velocity: {self.impact_velocity}
        impact_angle: {self.impact_angle}
        normal_velocity: {self.normal_velocity}
        """
