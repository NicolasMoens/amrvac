import numpy as np
import matplotlib.pyplot as plt
import os

c0 = 1.2727881099940257E+053   
c1 = 1.8869620128175814E+040

a = 6642027851021.7812        
b = 6642027851021.8145

N = 1000

def f(x):
    return x**4.0 + c1*x - c0

e = np.linspace(0,abs(c0/c1),N)

plt.plot(e, f(e))
plt.plot(a*np.ones(N), f(e))
plt.plot(b*np.ones(N), f(e))

plt.show() 