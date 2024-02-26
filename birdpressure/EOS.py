from birdpressure.Bird import Bird

import sympy as sy
import numpy as np


class EOS():
    """ Equation of state function to determine c_r

        Parameters
        ----------
        bird: class
        name: str
            name of eos

        Attributes
        ----------
        bird: class
        name: str
            name of eos
        name_lst: list
            list of available eos names
        rho: sympy.core.symbol.Symbol
            density sympu symbol
        P_0: int
            reference pressure Murnaghan [Pa]
        B: int
             Murnaghan material constant [Pa]
        gamma: float
            Murnaghan material constant [-]
        P_func: sympy.core.mul.Mul
            sympy function of chosen eos
        K: int
            linear eos bulk modulus [Pa]
        mu: sympy.core.mul.Mul
            change of density during impact
        k : int
            polynomial eos experimental constant
    """

    def __init__(self, bird: Bird, name: str):


        self.bird = bird
        self.name = name
        self.name_lst = ['Murnaghan', 'Linear_EOS', 'Polynomial_EOS']

        self.rho = sy.symbols('rho', real = True, positive=True)

        if self.name == "Murnaghan":
            self.P_0 = 0
            self.B = 128e6
            self.gamma = 7.98


            self.P_func = self.P_0 + self.B * ((self.rho / self.bird.get_density()) ** self.gamma - 1)


        if self.name == 'Linear_EOS':
            self.K = 2200e6
            self.P_func = self.K * (self.rho / self.bird.get_density()-1)

        if self.name == 'Polynomial_EOS':
            self.mu = self.rho/self.bird.get_density()-1
            self.k = 2

            c_1 = self.bird.get_density() * self.bird.get_sound_speed()**2
            c_2 = (2*self.k-1)*c_1
            c_3 = (self.k-1)*(3*self.k-1)*c_1

            self.P_func = c_1 * self.mu + c_2 * self.mu**2 + c_3 * self.mu**3

        if self.name == 'Wilbeck':
            pass




    def get_c_r(self,peak_pressure):
        """Compute release wave velocity by determining eos slope at peak pressure

         Parameters
         ----------
         peak_pressure: float
            hugonoit pressure [Pa]

        Returns
        -------
        c_r: float
            release wave velocity [m/s]
         """
        rho_peak = self.get_density(peak_pressure)

        c_r_squared = sy.diff(self.P_func)
        c_r = np.sqrt(float(c_r_squared.subs(self.rho,rho_peak)))
        return c_r

    def get_density(self,P):
        """Compute density from pressure"""
        eq = sy.Eq(self.P_func,P)

        return sy.nsolve(eq,self.rho, self.bird.get_density())

    def get_pressure(self,rho):
        """Compute pressure from density"""
        return self.P_func.subs(self.rho,rho)



    def eos_list(self):
        """Print string of eos options"""
        print('Current options of EOS:')
        for eos in self.name_lst:
            print(eos)