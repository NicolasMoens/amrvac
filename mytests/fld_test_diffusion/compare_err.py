import numpy as np
import matplotlib.pyplot as plt
import os

f1 = os.getcwd() + '/3x3error_out1d-5'
f2 = os.getcwd() + '/3x3error_out1d-4'
f3 = os.getcwd() + '/3x3error_out1d-3'
g1 = os.getcwd() + '/6x6error_out1d-5'
g2 = os.getcwd() + '/6x6error_out1d-4'
g3 = os.getcwd() + '/6x6error_out1d-3'

it_f1, err_f1 = np.loadtxt(f1, unpack='true')
it_f2, err_f2 = np.loadtxt(f2, unpack='true')
it_f3, err_f3 = np.loadtxt(f3, unpack='true')
it_g1, err_g1 = np.loadtxt(g1, unpack='true')
it_g2, err_g2 = np.loadtxt(g2, unpack='true')
it_g3, err_g3 = np.loadtxt(g3, unpack='true')

plt.plot(it_f1, err_f1, 'kx', label= '$30 \\times 30$ ; $Err = 10^{-5}$')
plt.plot(it_f2, err_f2, 'bx', label= '$30 \\times 30$ ; $Err = 10^{-4}$')
plt.plot(it_f3, err_f3, 'rx', label= '$30 \\times 30$ ; $Err = 10^{-3}$')

plt.plot(it_g1, err_g1, 'k+', label= '$60 \\times 60$ ; $Err = 10^{-5}$')
plt.plot(it_g2, err_g2, 'b+', label= '$60 \\times 60$ ; $Err = 10^{-4}$')
plt.plot(it_g3, err_g3, 'r+', label= '$60 \\times 60$ ; $Err = 10^{-3}$')

plt.title('ADI: Comparisson of residuals')
plt.xlabel('iteration', fontsize = '20')
plt.ylabel('mean residual', fontsize = '20')
plt.legend(loc = 'best')
plt.show()
