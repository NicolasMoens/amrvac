import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import os

nx = 50
ny = 100

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

f1 = os.getcwd() + '/const_2.okc'
f2 = os.getcwd() + '/const_6.okc'
f3 = os.getcwd() + '/const_10.okc'
f4 = os.getcwd() + '/kram_2.okc'
f5 = os.getcwd() + '/kram_6.okc'
f6 = os.getcwd() + '/kram_10.okc'
f7 = os.getcwd() + '/bump_2.okc'
f8 = os.getcwd() + '/bump_6.okc'
f9 = os.getcwd() + '/bump_10.okc'


#axi = (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)
fi = [f1, f2, f3, f4, f5, f6, f7, f8, f9]

f, ax = plt.subplots(3,3)
f.subplots_adjust(hspace =.2, wspace=.001)
ax = ax.ravel()

plt.setp(ax, xticks=[-0.75, 0.75])

mins = np.ones(9)*0.01
maxs = np.ones(9)*1

for i in range(len(fi)):
    x,y,z = get_val(fi[i])
    c = ax[i].pcolormesh(x, y, z, cmap = 'RdBu',norm=LogNorm( vmin = mins[i],vmax = maxs[i]))
    if (i == 2) or (i==5) or (i == 8):
        cbar = f.colorbar(c, ax=ax[i], ticks=[0.01, 0.1, 1])
    ax[i].set_aspect('1.')
    #ax[i].axis([-1,1,0,4])
    ax[i].axis([np.array(x).min(), np.array(x).max(), np.array(y).min(), np.array(y).max()])
    ax[i].set_xlim([np.array(x).min(), np.array(x).max()])
    ax[i].set_ylim([np.array(y).min(), np.array(y).max()])



    #ax[i].set_axis_off()

ax[1].set_title('Isothermal atmospheres, different $\kappa$ laws', fontsize=20)

ax[0].set_ylabel('Constant', fontsize=15)
ax[3].set_ylabel('Kramers', fontsize=15)
ax[6].set_ylabel('Bump', fontsize=15)

ax[6].set_xlabel('$t = 2 \\tau_{hd}$', fontsize=15)
ax[7].set_xlabel('$t = 6 \\tau_{hd}$', fontsize=15)
ax[8].set_xlabel('$t = 10 \\tau_{hd}$', fontsize=15)

# plt.tight_layout()

plt.savefig('density_kappa_grid')
plt.show()
