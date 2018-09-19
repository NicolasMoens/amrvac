import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import os

nx = 20
ny = 40

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

f1 = os.getcwd() + '/iso_0.okc'
f2 = os.getcwd() + '/iso_2.okc'
f3 = os.getcwd() + '/iso_4.okc'
f4 = os.getcwd() + '/iso_6.okc'
f5 = os.getcwd() + '/iso_8.okc'
f6 = os.getcwd() + '/iso_10.okc'

# g1 = os.getcwd() + '/iso_0.okc'
# g2 = os.getcwd() + '/iso_2.okc'
# g3 = os.getcwd() + '/iso_4.okc'
# g4 = os.getcwd() + '/iso_6.okc'
# g5 = os.getcwd() + '/iso_8.okc'
# g6 = os.getcwd() + '/iso_10.okc'

g1 = os.getcwd() + '/non_0.okc'
g2 = os.getcwd() + '/non_2.okc'
g3 = os.getcwd() + '/non_4.okc'
g4 = os.getcwd() + '/non_6.okc'
g5 = os.getcwd() + '/non_8.okc'
g6 = os.getcwd() + '/non_10.okc'


#axi = (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)
fi = [f1, f2, f3, f4, f5, f6, g1, g2, g3, g4, g5, g6]

f, ax = plt.subplots(2,6, sharex = True, sharey = True)
f.subplots_adjust(hspace =.1, wspace=.0)
ax = ax.ravel()

mins = np.ones(12)*0.01
maxs = np.ones(12)*1.0

for i in range(len(fi)):
    x,y,z = get_val(fi[i])
    c = ax[i].pcolor(x, y, z, cmap = 'RdBu',norm=LogNorm( vmin = mins[i],vmax = maxs[i]))
    if (i == 5) or (i==11):
        cbar = f.colorbar(c, ax=ax[i], ticks=[0.01, 0.1, 1])
    ax[i].set_aspect('1.')
    ax[i].axis([-1,1,0,4])
    ax[i].set_axis_off()

ax[3].set_title('Isothermal and non-isothermal atmospheres', fontsize=20)

ax[0].set_ylabel('Isothermal', fontsize=15)
ax[6].set_ylabel('Non-Iso', fontsize=15)

ax[6].set_xlabel('$t = 0t_0$', fontsize=15)
ax[7].set_xlabel('$t = 2t_0$', fontsize=15)
ax[8].set_xlabel('$t = 4t_0$', fontsize=15)
ax[9].set_xlabel('$t = 6t_0$', fontsize=15)
ax[10].set_xlabel('$t = 8t_0$', fontsize=15)
ax[11].set_xlabel('$t = 10t_0$', fontsize=15)


# plt.tight_layout()

plt.savefig('density_kappa_grid')
plt.show()
