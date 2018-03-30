import numpy as np
import matplotlib.pyplot as plt

M_sun = 1.99e33 #cgs
L_sun = 3.9e33 #cgs
R_sun = 6.96e10

mu = 0.6 #cgs
m_p = 1.67e-24 #cgs
k = 1.38e-16 #cgs
G = 6.67e-8
kappa = 0.34
c_light = 2.99e10
gamma = 5.0/3.0

######################################
M_star = 50*M_sun #cgs
L_star = 1000*L_sun #cgs
T_star = 50000 #cgs
R_star = 1000*R_sun 
rho_0 = 1e2
######################################

F = L_star/(4*np.pi*R_star**2)

c_sound = np.sqrt(k*T_star/(mu*m_p))

E_0 = c_sound*T_star**4
p_0 = rho_0*c_sound**2


def IntegrateInitCond(E_0, p_0, R_max):
    E = [E_0]
    p = [p_0]
    e = [p_0*1.0/(gamma - 1.0)]
    dx = 0.01*R_sun
    r = np.linspace(R_star,R_star+R_max,R_max/dx)
    for x in r[1:]:
        dp = (rho_0*kappa/c_light)*F - rho_0*(G*M_star/(R_star**2))
        p.append(p[-1] + dx*dp)
        e.append(p[-1]*1.0/(gamma - 1.0))
        dE = -3*(rho_0*kappa/c_light)*F
        E.append(E[-1] + dx*dE)
    return r,E,e

r,E,e = IntegrateInitCond(E_0,p_0, 1.0*R_sun)

print len(r)
print len(E)
print len(e)

plt.figure()
#plt.plot(r,E,label='Radiative Energy')
plt.plot(r,e,label='Gas Energy')
plt.legend()
plt.show()
        
