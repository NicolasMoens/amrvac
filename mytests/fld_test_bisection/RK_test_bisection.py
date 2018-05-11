import numpy as np
import matplotlib.pyplot as plt
import os

E = 1000000000000.0000
c = 29979245800.000000 
rho = 9.9999999999999995e-008
kappa = 0.40000000000000002 
sigma = 1.6767990825413562e-037
gamma = 1.6666669999999999 

e0 = rho/(gamma - 1.)*(c*E/(4.*sigma))**(1./4.)

print e0

def f(e_i):
    a = c*rho*kappa*E
    b = 4.*rho*kappa*sigma*((gamma-1)/rho)**4.0
    
    return ( a - b*e_i**4.0 )
    

def GetEnergyProfile(frac,N):
        
    t_0 = 0
    t_f = 1.e-7
    dt = (t_f-t_0)/N
        
    t = np.linspace(t_0,t_f,N)

    e_f = frac*e0

    e = []
    
    for i in t:
        k1 = f(e_f)
        k2 = f(e_f+dt/2.0*k1)
        k3 = f(e_f+dt/2.0*k2)
        k4 = f(e_f+dt*k3)

        e_f = e_f + dt*(k1 + 2*k2 + 2*k3 + k4)/6.0
        e.append(e_f/e0)

    plt.loglog(t,e,'k--', label='explicit $e_{gas}/e_{eq}$ = ' + str(frac))

    
    print "RK model done"
    return


t_0 = 0
t_f = 1.e-7
x = np.linspace(t_0,t_f,1000000)
plt.loglog(x,np.ones(len(x))*e0/e0,'b', label = 'radiative equilibrium')

for frac in [1.e-10,1.e-4,1.e2]: 
                                                                                                                                                                                                                                                                GetEnergyProfile(frac,1000000)

                                                                                                                                                                                                                                                                f1 = os.getcwd() + '/energy_out1'
f2 = os.getcwd() + '/energy_out2'
f3 = os.getcwd() + '/energy_out3'
#f4 = os.getcwd() + '/energy_out4'
#f5 = os.getcwd() + '/energy_out5'

for f in [f1,f2,f3]:
    it,t,e, r_e = np.loadtxt(f, unpack='true')
    plt.loglog(t,e/e0, 'r+', label = 'implicit')
    print "implicit model done"

plt.legend(loc = 'best')
plt.ylabel('$\log(e_{gas}/e_{eq})$')
plt.xlabel('$\log(t[s])$')
plt.show()
