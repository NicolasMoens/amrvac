import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import os

nx = 40
ny = 80

v0 = 2491838.8019151785
T0 = 45788.310074951703
e0 = 556097.85659898631
t0 = 1317.9222983369184
rho0 = 8.9559432451549764e-008
l0 = 3284049920.9051652

kappa = 0.34*v0*t0*rho0
sigma = 5.67e-5*T0**4/(e0*v0)
hd_gamma = 1.6666667
speedoflight = 2.99e10/v0
mp = 1.672621777e-24/rho0*l0**3
kb = 1.3806488e-16/rho0*l0**3*T0

dt = 5.601E-04

def get_val(f):
    X, Y, Z, R, EG, ER, Temp = np.loadtxt(f, unpack='true', skiprows = 15)
    x = []
    y = []
    rho = []
    e_g = []
    e_r = []
    tmp = []
    for i in range(ny):
        x.append(X[i*nx:(i+1)*nx])
        y.append(Y[i*nx:(i+1)*nx])
        rho.append(R[i*nx:(i+1)*nx])
        e_g.append(EG[i*nx:(i+1)*nx])
        e_r.append(ER[i*nx:(i+1)*nx])
        tmp.append(ER[i*nx:(i+1)*nx])
    return x,y,rho,e_g,e_r,tmp

f1 = os.getcwd() + '/noniso_const.okc'
f1 = os.getcwd() + '/noniso_test.okc'

x,y,rho,e_g,E_r,tmp = get_val(f1)

heating = np.ones([ny,nx])
cooling = np.ones([ny,nx])
exch = np.ones([ny,nx])

for i in range(nx):
    for j in range(ny):

        # rho[j][i] = rho[j][i]*rho0
        # e_g[j][i] = e_g[j][i]*e0
        # E_r[j][i] = E_r[j][i]*e0
        # tmp[j][i] = tmp[j][i]*T0

        heating[j][i] = 4.*kappa*rho[j][i]*sigma*((hd_gamma-1.)*(e_g[j][i])/(rho[j][i]))**4.
        cooling[j][i] = speedoflight*kappa*rho[j][i]*E_r[j][i]

        # heating[j][i] = 4.*sigma*((hd_gamma-1.)*(e_g[j][i])/(rho[j][i]))**4.
        # cooling[j][i] = speedoflight*E_r[j][i]

        exch[j][i] = (heating[j][i] - cooling[j][i])*dt

plt.plot(heating[:][5])
plt.plot(cooling[:][5])
plt.show()



f, ax = plt.subplots(2,2)
f.subplots_adjust(hspace =.6, wspace=.001)
ax = ax.ravel()

plt.setp(ax, xticks=[-0.75, 0.75])
a = ax[2].pcolormesh(x, y, rho, cmap = 'RdBu', vmin = np.array(rho).min(),vmax = np.array(rho).max())
# b = ax[0].pcolormesh(x, y, e_g, cmap = 'RdBu', vmin = np.array(e_g).min(),vmax = np.array(e_g).max())
b = ax[0].pcolormesh(x, y, tmp, cmap = 'RdBu', vmin = np.array(tmp).min(),vmax = np.array(tmp).max())
c = ax[1].pcolormesh(x, y, E_r, cmap = 'RdBu', vmin = np.array(E_r).min(),vmax = np.array(E_r).max())
d = ax[3].pcolormesh(x, y, exch, cmap = 'RdBu', vmin = np.array(exch).min(),vmax = np.array(exch).max())

abar = f.colorbar(a, ax=ax[2]) #, ticks=[np.array(rho).min(), np.array(rho).max()])
bbar = f.colorbar(b, ax=ax[0]) #, ticks=[np.array(e_g).min(), np.array(e_g).max()])
cbar = f.colorbar(c, ax=ax[1]) #, ticks=[np.array(E_r).min(), np.array(E_r).max()])
dbar = f.colorbar(d, ax=ax[3]) #, ticks=[np.array(exch).min(), np.array(exch).max()])

for i in range(4):
        ax[i].set_aspect('1.')
        ax[i].axis([np.array(x).min(), np.array(x).max(), np.array(y).min(), np.array(y).max()])
        ax[i].set_xlim([np.array(x).min(), np.array(x).max()])
        ax[i].set_ylim([np.array(y).min(), np.array(y).max()])

print 'heating'
print heating[:][0]
print 'cooling'
print cooling[:][0]

for i in range(nx):
    print i,heating[0][i], cooling[0][i], heating[0][i]/cooling[0][i], heating[0][i]-cooling[0][i], exch[0][i]


ax[2].set_title('$\\rho/\\rho_0$', fontsize = 15)
ax[0].set_title('$T/T_0$', fontsize = 15)
#ax[0].set_title('$e/e_0$', fontsize = 15)
ax[1].set_title('$E/e_0$', fontsize = 15)
ax[3].set_title('$c \kappa \\rho E - 4 \\pi \kappa \\rho B(T)$ \n $e_0/dt$', fontsize = 15)

plt.savefig('E_exchange')
plt.show()
