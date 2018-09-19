import numpy as np
import matplotlib.pyplot as plt
import os

f0 = os.getcwd() + '/energy_out0'
f1 = os.getcwd() + '/energy_out1'
f2 = os.getcwd() + '/energy_out2'
f3 = os.getcwd() + '/energy_out3'

F = [f0]

for f in F:
    i = np.loadtxt(f,usecols=[0])
    print 'file read'
    print i
    t = np.loadtxt(f,usecols=[1])
    print 'file read'
    print t
    e = np.loadtxt(f,usecols=[2])
    print 'file read'
    print e

    plt.loglog(t,e,'+',label='gas-energy')
    
    print 'plotted'
    
#plt.loglog(t,1.14950000e-05*np.ones(len(t)), 'r--', label='$4 \pi \kappa B = c \kappa E$')
plt.title('Convergence to radiative Equilibrium', fontsize=20)
#plt.xlabel('log(t/t0)', fontsize=20)
#plt.ylabel('log(e/e0)', fontsize=20)
#plt.ylim([1e-7,1e-4])
#plt.legend(loc='lower right', fontsize=20)
plt.show()
