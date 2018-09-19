import numpy as np
import matplotlib.pyplot as plt
import os

f1 = os.getcwd() + '/const_i.okc'
f2 = os.getcwd() + '/const_f.okc'

nx = 20
ny = 40

y = np.linspace(0,4,ny)

def get_mean(f):
    G, K = np.loadtxt(f, unpack='true', skiprows = 5)
    g_mean = []
    g_var = []
    k_mean = []
    k_var = []
    for i in range(ny):
        g_m = np.mean(G[i*nx:(i+1)*nx])
        g_v = np.std(G[i*nx:(i+1)*nx])
        g_mean.append(g_m)
        g_var.append(g_v)
        k_m = np.mean(K[i*nx:(i+1)*nx])
        k_v = np.std(K[i*nx:(i+1)*nx])
        k_mean.append(k_m)
        k_var.append(k_v)
    return g_mean, g_var, k_mean, k_var

rhoI, rhoI_v, pI, pI_v = get_mean(f1)
rhoF, rhoF_v, pF, pF_v = get_mean(f2)

print len(y)
print len(rhoI)
print len(rhoI_v)
print  rhoI_v

rhoI_p = np.array(rhoI) + np.array(rhoI_v)
rhoI_m = np.array(rhoI) - np.array(rhoI_v)
rhoF_p = np.array(rhoF) + np.array(rhoF_v)
rhoF_m = np.array(rhoF) - np.array(rhoF_v)

pI_p = np.array(pI) + np.array(pI_v)
pI_m = np.array(pI) - np.array(pI_v)
pF_p = np.array(pF) + np.array(pF_v)
pF_m = np.array(pF) - np.array(pF_v)

f, (ax1, ax2) = plt.subplots(2,1, sharex = True)
ax1.set_title('$\\rho$ and $R_{rad}$, initial and final state', fontsize = 20)
ax1.plot(y,rhoI,'r-', label = 'Initial')
ax1.plot(y,rhoI_p, 'r--')
ax1.plot(y,rhoI_m, 'r--')
ax1.plot(y,rhoF,'b-', label = 'Final')
ax1.plot(y,rhoF_p, 'b--')
ax1.plot(y,rhoF_m, 'b--')
ax1.set_ylabel('$\\rho/ \\rho_0$', fontsize = 15)
ax1.legend(loc = 'best', fontsize = 15)

ax2.plot(y,pI,'r-', label = 'Initial')
ax2.plot(y,pI_p, 'r--')
ax2.plot(y,pI_m, 'r--')
ax2.plot(y,pF,'b-', label = 'Final')
ax2.plot(y,pF_p, 'b--')
ax2.plot(y,pF_m, 'b--')
ax2.set_xlabel('$y/H_{eff}$', fontsize = 15)
ax2.set_ylabel('$E_{rad}/e_0$', fontsize = 15)

f.savefig(subplot_density)


plt.figure(1)
plt.title('$\\rho$ in initial and final state', fontsize = 20)
plt.plot(y,rhoI,'r-', label = 'Initial')
plt.plot(y,rhoI_p, 'r--')
plt.plot(y,rhoI_m, 'r--')
plt.plot(y,rhoF,'b-', label = 'Final')
plt.plot(y,rhoF_p, 'b--')
plt.plot(y,rhoF_m, 'b--')
plt.ylabel('$\\rho/ \\rho_0$', fontsize = 15)
plt.xlabel('$y/H_{eff}$', fontsize = 15)
plt.legend(loc = 'best', fontsize = 15)
plt.savefig(density)

plt.show()
