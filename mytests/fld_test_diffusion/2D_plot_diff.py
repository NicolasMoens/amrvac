import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

nx = 30
ny = 30

def get_val(f):
    X, Y, Z, F = np.loadtxt(f, unpack='true', skiprows = 9)
    x = []
    y = []
    f = []
    for i in range(ny):
        x.append(X[i*nx:(i+1)*nx])
        y.append(Y[i*nx:(i+1)*nx])
        f.append(F[i*nx:(i+1)*nx])
    return x,y,f

f1 = os.getcwd() + '/num_0.okc'
f2 = os.getcwd() + '/num_5.okc'
f3 = os.getcwd() + '/num_10.okc'
f4 = os.getcwd() + '/theo_0.okc'
f5 = os.getcwd() + '/theo_5.okc'
f6 = os.getcwd() + '/theo_10.okc'
f7 = os.getcwd() + '/res_0.okc'
f8 = os.getcwd() + '/res_5.okc'
f9 = os.getcwd() + '/res_10.okc'


#axi = (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)
fi = [f1, f2, f3, f4, f5, f6, f7, f8, f9]

f, ax = plt.subplots(3,3, sharex = True, sharey = True)
f.subplots_adjust(hspace =.1, wspace=.001)
ax = ax.ravel()


mins = [1, 1, 1, 1, 1, 1, 0.00001, 0.00001, 0.00001]
maxs = [3, 3, 3, 3, 3, 3, 0.01, 0.01, 0.01]

for i in range(len(fi)):
    x,y,z = get_val(fi[i])
    c = ax[i].pcolor(x,y,z, cmap = 'RdBu',vmin = mins[i],vmax = maxs[i])
    ax[i].axis([-0.5,0.5,0,1])
    if (i == 2) or (i==5):
        cbar = f.colorbar(c, ax=ax[i], ticks=[1, 2, 3])
    if (i == 8):
        cbar = f.colorbar(c, ax=ax[i], ticks=[0.00001, 0.01])

ax[1].set_title('Radiative diffusion', fontsize=20)

ax[0].set_ylabel('Numerical', fontsize=15)
ax[3].set_ylabel('Theoretical', fontsize=15)
ax[6].set_ylabel('Residual', fontsize=15)

ax[6].set_xlabel('$t = 0t_0$')
ax[7].set_xlabel('$t = 0.5t_0$', fontsize=15)
ax[8].set_xlabel('$t = 1t_0$', fontsize=15)

# ax[2].set_colorbar()
# ax[5].set_colorbar()
# ax[8].set_colorbar()

plt.show()
