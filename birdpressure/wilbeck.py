import numpy as np
import sympy as sp

#Linear Hugoniot
##  in: c_0[array/float] = isentropic wave speed, k[float] = material const., u_0[array/float] = impact velociuty
## out: u_s[array/float] = shock velocity
def find_u_s(c_0, k, u_0):

    u_s = c_0 + k * u_0

    return u_s

#Finding speed of sound in shocked medium with Murnaghan EOS
##  in: c_0[array/float] = isentropic wave speed, k[float] = material const., u_0[array/float] = impact velocity
## out: c_r[array/float] =  speed of sound in shocked medium
def find_c_r(c_0,k ,u_0):
    u_s     = find_u_s(c_0, k, u_0)
    rho_0   = 950
    rho     = rho_0 * u_s / (u_s - u_0)
    B       = 128e6
    lam     = 7.89

    c_r     = np.sqrt(lam*B/rho_0**lam * rho**(lam-1)) +451.8 #added so u_s and c_r start at same point
    return c_r

#Determining Peak pressure Impact Duaration
##  in: c_r[array/float] =  speed of sound in shocked medium, a[array/float] = impactor radius
## out: t_b[array/float] = P_H duaration
def find_t_b(c_r,a):
    t_b = a/c_r
    return t_b

#Determine t_c
##  in: a[array/float] = impactor radius, c_r[array/float] =  speed of sound in shocked medium, u_s[array/float] = shock velocity, u_0[array/float] = impact velociuty
## out: t_c[array/float] =  time release wave needs to initiate steady state
def find_t_c(a,c_r,u_s,u_0):
    t_c = a / np.sqrt(c_r ** 2 - (u_s - u_0) ** 2)
    return t_c

# Determine Critical projectile aspect ratio L/D_c
##  in:  u_s[array/float] = shock velocity, c_r[array/float] = speed of sound in shocked medium, u_0[array/float] = impact velociuty
## out: LC[array/float] = critical bird aspect ratio
def find_LD_crit(u_s,c_r,u_0):

    if len(u_s)>1:
        LC      = np.zeros_like(u_s)
        LC[1::] = u_s[1::]/(2*np.sqrt(c_r[1::]**2-(u_s[1::]-u_0[1::])**2))
        LC[0]   = 10e5  # avoid devision by 0 as if u_0 = 0 than c_r = u_s
    else:

        if u_0 ==0:
            LC = 10e100
        else:
            LC = u_s/(2*np.sqrt(c_r**2-(u_s-u_0)**2))

    return LC


# Banks & Chandrasekhara radial pressure distribution (changing rho)
##  in: u_0[float] = impact velocity, r[float] = distance from center, P_start[float] = initial pressure,
#       rho_start[float]= initial density
## out: res[array(1,2)] = array with resulting density and pressure [rho,P]
def solve_Banks(u_0,r,P_start,rho_start):
    rho_0   = 950
    B       = 128e6
    gamma   = 7.98
    P_0     = 101325
    xi_1    = 0.5
    a       = 0.5

    rho = sp.Symbol('rho', real = True, positive = True)
    P   = sp.Symbol('P', real = True, positive = True)

    eq1 = sp.Eq(P, P_0 + B * ((rho / rho_0) ** gamma - 1))
    eq2 = sp.Eq(P, 1/2 * rho * u_0**2 * np.e**(-xi_1 * (r/a)**2))

    res = sp.nsolve((eq1,eq2),(rho,P),(rho_start,P_start))
    return res

# Banks & Chandrasekhara radial pressure distribution (const. rho)
##  in: rho[float] = density, u_0[float] = impact velocity,xi_1[float]= material const.,
#       r[array/float] = distance from center, a[float]= bird radius
## out: P[array/float] = Pressure
def find_banks(rho,u_0,xi_1,r,a):
    P = 1/2 * rho * u_0**2 * np.e**(-xi_1 * (r/a)**2)
    return P

# Leach & Walker radial pressure distribution (changing rho)
##  in: u_0[float] = impact velocity, r[float] = distance from center, P_start[float] = initial pressure,
#       rho_start[float]= initial density
## out: res[array(1,2)] = array with resulting density and pressure [rho,P]
def solve_Leach(u_0,r,P_start,rho_start):
    rho_0   = 950
    B       = 128e6
    gamma   = 7.98
    P_0     = 101325
    xi_2    = 2.58
    a       = 0.5

    rho = sp.Symbol('rho', real = True, positive = True)
    P   = sp.Symbol('P', real = True)

    eq1 = sp.Eq(P, P_0 + B * ((rho / rho_0) ** gamma - 1))
    eq2 = sp.Eq(P, 1/2 * rho * u_0**2 * (1 - 3 * (r/(xi_2*a))**2 + (r/(xi_2*a))**3))

    res = sp.nsolve((eq1,eq2),(rho,P),(rho_start,P_start))
    return res

# # Leach & Walker radial pressure distribution (const. rho)
##  in: rho[float] = density, u_0[float] = impact velocity,xi_2[float]= material const.,
#       r[array/float] = distance from center, a[float]= bird radius
## out: P[array/float] = Pressure
def find_leach(rho, u_0, xi_2, r, a):
    P = 1 / 2 * rho * u_0 ** 2 * (1 - 3 * (r/(xi_2*a)**2 +(r/(xi_2*a)**2 )**3))
    return P

# Murnaghan EOS
##  in: rho[array/float] = density
## out: P[array/float]  = pressure
def eos(rho):
    rho_0   = 950
    P_0     = 101325
    B       = 128e6
    gamma   = 7.98
    P       = P_0 + B * ( (rho/rho_0)**gamma - 1)
    return P
# find position of nearest value
##  in:  array[array] = numpy array, value[float] = wanted value
## out: idx[int] = index of closest value in array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx