import numpy as np
import matplotlib.pyplot as plt

c0 = -4e-4
c1 = -2e-2

xmax = max(abs(c0/c1),abs(c0)**(1.0/4.0))

x = np.linspace(-0.2,1.2*xmax,1000)

def poly(x,c0,c1):
    return x**4.0 + c1*x - c0 

def bisect(c0_l,c1_l):
    lb = 0
    ub = xmax
    cb = ub
    
    while abs(lb - ub) > 1e-6:
        
        print poly(lb,c0_l,c1_l), poly(cb,c0_l,c1_l), poly(ub,c0_l,c1_l)
        
        if poly(lb,c0_l,c1_l)*poly(cb,c0_l,c1_l) < 0:
            ub = cb 
        elif poly(ub,c0_l,c1_l)*poly(cb,c0_l,c1_l) < 0:
            lb = cb
        else:
            print "ERROR"
            exit()
        cb = (ub + lb)/2.0
    print cb
    return cb

n0 = 1000
n1 = 1000

c0_arr = np.linspace(1.9,10,n0)
c1_arr = np.linspace(1.2,10,n1)

C0,C1 = np.meshgrid(c0_arr,c1_arr)
root = np.zeros((n0,n1))

for i in range(n0):
    for j in range(n1):
        root[i][j] = bisect(C0[i,j], C1[i,j])
        
