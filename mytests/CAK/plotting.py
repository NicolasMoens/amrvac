import numpy as np
import matplotlib.pyplot as plt
import os


f = os.getcwd() + '/test_out'

R = 1.39e12

i = np.loadtxt(f,usecols=[0])
r = np.loadtxt(f,usecols=[1]) #*R
v = np.loadtxt(f,usecols=[2]) #*2.34e6
rho = np.loadtxt(f,usecols=[3])
g_cak = np.loadtxt(f,skiprows=2,usecols=[4])
g_sc_grav = np.loadtxt(f,skiprows=2,usecols=[5])
dvdr = np.loadtxt(f,usecols=[6])
fd = np.loadtxt(f,usecols=[7])


g = os.getcwd() + '/test_init'
i_i = np.loadtxt(g,usecols=[0])
r_i = np.loadtxt(g,usecols=[1]) #*R
v_i = np.loadtxt(g,usecols=[2]) #*2.34e6
rho_i = np.loadtxt(g,usecols=[3])

M_dot = 4*np.pi*np.power(r,2)*rho*v
M_dot_i =4*np.pi*np.power(r_i,2)*rho_i*v_i

M_dot = M_dot #*(365*24*60*60)/(2.0e33)
M_dot_i = M_dot_i #*(365*24*60*60)/(2.0e33)


print '-----------------------------------------------'
print 'RADIUS', r[-1]
print 'velocity', v[0:20]
print 'rho',rho[0:20]


kappa_e = 0.34
L =8e5*3.9e33
q = 2000.
alpha = 0.67
gamma_e = 0.34
G = 6.6e-8
M = 50*1.9e33
escape_speed = (2.*G*M/R)**0.5




#M_dot_ = 4*1e-6/(365*24*60*60)*(1.9e33)
g_e = kappa_e*L/(4*np.pi*np.power(r,2)*2.99e10)

gamma_e = g_e*np.power(r,2)/(G*M)
v_inf = (1-gamma_e)**0.5*escape_speed*(alpha/(1-alpha))**0.5

g_th_cak =  g_e*q/(1-alpha)*((4*np.pi*R*v_inf)/(q*M_dot*kappa_e))**alpha

t1 = kappa_e*rho*2.99e10/dvdr
g_th_1 = q*gamma_e/((1-alpha)*(q*t1)**alpha)*fd

t2 = q*rho*kappa_e*2.99e10/(R*v_inf**2/(np.power(r,2)*v*23e5))
g_th_2 = (q*g_e)/(t2**alpha)*fd

print 'v_inf', v_inf
print 'g_e', g_e[3:9]
print 'g_th_cak', g_th_cak[3:9]
print 'g_cak', g_cak[3:9]

plt.figure(1)
plt.title('Velocity')
plt.plot(r,v,label= 'v')
plt.plot(r_i,v_i,'.', label= 'v_init')
plt.legend()

plt.figure(2)
plt.title('Density')
plt.semilogy(r,rho,label = 'rho')
plt.semilogy(r_i,rho_i,label = 'rho_init')
plt.legend()

plt.figure(3)
plt.title('Massloss')
plt.plot(r,M_dot,label = 'Massloss')
plt.plot(r_i,M_dot_i,label = 'Massloss_init')

plt.figure(5)
plt.title('finite disc correction')
plt.plot(r[:-2],fd[:-2])

plt.figure(4)
plt.title('accelerations')
plt.plot(r[2:],g_cak,label = 'cak')
plt.plot(r[2:],g_sc_grav,label = 'grav')
#plt.plot(r,g_th_cak,label = 'g_th_cak')
plt.plot(r,g_th_1,label = 'g_th_1')
plt.plot(r,g_th_2,label = 'g_th_2')
plt.plot(r,g_e,label = 'g_e')
plt.legend()

plt.legend()
plt.show()
