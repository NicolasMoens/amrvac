import numpy as np
import matplotlib.pyplot as plt
import os

f1 = os.getcwd() + '/const.okc'
f2 = os.getcwd() + '/kramers.okc'
f3 = os.getcwd() + '/bump.okc'

unit_kappa = 0.34/100

nx = 20
ny = 40

y = np.linspace(0,4,ny)


def get_mean(f):
    K, G = np.loadtxt(f, unpack='true', skiprows = 5)
    g_mean = []
    k_mean = []
    for i in range(ny):
        g_m = np.mean(G[i*nx:(i+1)*nx])
        g_mean.append(g_m)
        k_m = np.mean(K[i*nx:(i+1)*nx])
        k_mean.append(k_m*unit_kappa)
    return g_mean, k_mean

g_const, k_const = get_mean(f1)
g_kramers, k_kramers = get_mean(f2)
g_bump, k_bump = get_mean(f3)


f, (ax1, ax2) = plt.subplots(2,1, sharex = True)
ax1.set_title('$\Gamma$ and $\\kappa$ for different models', fontsize = 20)
ax1.plot(y,g_const, label = 'Constant $\\kappa$')
ax1.plot(y,g_kramers, label = 'Kramers $\\kappa$ law')
ax1.plot(y,g_bump, label = '$\\kappa$ bump')
ax1.set_ylabel('$\\Gamma$', fontsize = 15)
ax1.legend(loc = 'best', fontsize = 15)

ax2.plot(y,k_const, label = 'Constant $\\kappa$')
ax2.plot(y,k_kramers, label = 'Kramers $\\kappa$ law')
ax2.plot(y,k_bump, label = '$\\kappa$ bump')
ax2.set_xlabel('$y/H_{eff}$', fontsize = 15)
ax2.set_ylabel('$\\kappa [\\frac{cm^2}{g}]$', fontsize = 15)


plt.show()
