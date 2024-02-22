from BirdPressure.Bird import Bird

import sympy as sy
import numpy as np


class EOS():

    def __init__(self, bird: Bird, name: str):
        """
        Equation of state Class uses sympy to use different equations
        :param bird: Bird Class
        :param name: Name EOS you want to use ['Murnaghan', 'Linear_EOS', 'Polynomial_EOS']
        :param peak_pressur: Set pressure around which c_r is determined
        """


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
        rho = self.get_density(peak_pressure)

        c_r_squared = sy.diff(self.P_func)
        c_r = np.sqrt(float(c_r_squared.subs(self.rho,rho)))
        return c_r

    def get_density(self,P):
        eq = sy.Eq(self.P_func,P)

        return sy.nsolve(eq,self.rho, self.bird.get_density())

    def get_pressure(self,rho):
        return self.P_func.subs(self.rho,rho)



    def eos_list(self):
        print('Current options of EOS:')
        for eos in self.name_lst:
            print(eos)