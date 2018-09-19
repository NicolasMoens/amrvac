import numpy as np
import matplotlib.pyplot as plt

M_sun = 1.989e33
L_sun = 3.827e33
R_sun = 6.96e10

kb = 1.38e-16
sigma = 5.67e-5
mu = 0.6
mp = 1.67e-24
G = 6.67e-8
kappa = 0.34
c_light = 2.99e10
hd_gamma = 5.0/3.0


M_star = 150*M_sun
L_star = (M_star/M_sun)**3.e0 * L_sun
R_star = 30*R_sun
T_star = (L_star/(4*np.pi*sigma*R_star**2.0))**0.25

tau_0 = 100.0
T_bound = (1.0+(3.0/4.0*tau_0)**0.25)*T_star

################################################################################
####### ISOTHERMAL ATMOSPHERE
################################################################################

c_sound = np.sqrt(kb*T_star/(mu*mp))
g_grav = G*M_star/R_star**2.0

Flux = L_star/(4*np.pi*R_star**2.0)
Gamma = kappa*Flux/(c_light*g_grav)
g_eff = g_grav*(1.0-Gamma)
H_eff = c_sound**2.0/g_eff

p_bound = g_eff*tau_0/kappa
rho_bound = p_bound/c_sound**2.0
################################################################################
y = np.linspace(0,15*H_eff,1000)

p_iso = p_bound*np.exp(-y/H_eff)

rho_iso = 1.0/c_sound**2.0 * p_iso
e_iso = 1.0/(hd_gamma-1.0) * p_iso
T_iso = T_bound*np.ones(len(y))
E_iso = 3.0*Gamma/(1.0-Gamma) * p_iso

################################################################################
####### NON-ISOTHERMAL ATMOSPHERE
################################################################################

a = Flux*kappa/(4.0/3.0*sigma*g_eff)
b = T_bound**4.0-a*p_bound
c = -g_eff*mp*mu/kb

print "a", a, "b", b, "c", c

p_non = p_iso

for i in range(1,len(y)):
    k1 = (y[i]-y[i-1])*c*p_non[i-1]/(a*p_non[i-1] + b)**0.25
    k2 = (y[i]-y[i-1])*c*(p_non[i-1]+0.5*k1)/(a*(p_non[i-1]+0.5*k1) + b)**0.25
    k3 = (y[i]-y[i-1])*c*(p_non[i-1]+0.5*k2)/(a*(p_non[i-1]+0.5*k2) + b)**0.25
    k4 = (y[i]-y[i-1])*c*(p_non[i-1]+k3)/(a*(p_non[i-1]+k3) + b)**0.25

    p_non[i] = p_non[i-1] + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)

rho_non = mp*mu/kb * p_non/(a*p_non+b)**0.25
e_non = p_non/(hd_gamma-1.0)
T_non = (a*p_non+b)**(1.0/4.0)
E_non = 4.0*sigma/c_light*T_non**4.0

div_E = np.ones(len(E_non))
for i in range(len(div_E)-2):
    div_E[i+1] = (E_non[i+2] - E_non[i])/(y[i+2]-y[i])
# div_E[0] = (E_non[1] - E_non[0])/(y[1]-y[0])
# div_E[-1] = (E_non[-1] - E_non[-2])/(y[-1]-y[-2])
div_E[0] = div_E[1]
div_E[-1] = div_E[-2]

for i in range(len(div_E)):
    print i, div_E[i], rho_non[i], div_E[i]/rho_non[i],-c_light/(3.0*rho_non[i]*kappa)*div_E[i]

F_non = -c_light/(3.0*rho_non*kappa)*div_E

################################################################################
####### CHECK
################################################################################

#Redefine tau_0
tau_0 = 0
for i in range(len(y)-1):
    t1 = -(y[i+1]-y[i])*kappa*np.interp(y[i],y,rho_non)
    t2 = -(y[i+1]-y[i])*kappa*np.interp(y[i]+t1/2,y,rho_non)
    t3 = -(y[i+1]-y[i])*kappa*np.interp(y[i]+t2/2,y,rho_non)
    t4 = -(y[i+1]-y[i])*kappa*np.interp(y[i]+t3,y,rho_non)

    tau_0 = tau_0 - 1.0/6.0*(t1 + 2*t2 + 2*t3 + t4)

tau = tau_0*np.ones(len(y))

for i in range(len(y)-1):
    t1 = -(y[i+1]-y[i])*kappa*np.interp(y[i],y,rho_non)
    t2 = -(y[i+1]-y[i])*kappa*np.interp(y[i]+t1/2,y,rho_non)
    t3 = -(y[i+1]-y[i])*kappa*np.interp(y[i]+t2/2,y,rho_non)
    t4 = -(y[i+1]-y[i])*kappa*np.interp(y[i]+t3,y,rho_non)

    tau[i+1] = tau[i] + 1.0/6.0*(t1 + 2*t2 + 2*t3 + t4)
    print 1.0/6.0*(t1 + 2*t2 + 2*t3 + t4)


T_check = (3.0*Flux/(4.0*sigma)*(tau - tau_0) + T_bound**4.0)**0.25
T_check2 = (mp*mu)/kb*p_non/rho_non

E_check = 4.0*sigma/c*(a*p_non+b)

dEdt = np.ones(len(E_non))
for i in range(len(dEdt)-2):
    dEdt[i+1] = (E_non[i+2] - E_non[i])/(tau[i+2]-tau[i])

F_check = c_light/3.0*dEdt

################################################################################
####### Read from initial_cond
################################################################################

y_amr, rho_amr, e_amr, r_e_amr, T_amr, tau_amr = np.loadtxt('initial_cond', unpack = 'true')

################################################################################
####### Plotting
################################################################################

# for i in range(len(y)):
#     print i, y[i]/H_eff, tau[i], T_non[i], T_check[i], T_bound

################################################################################

# plt.figure(1)
# plt.title('Density $\\rho$')
# plt.plot(y_amr,rho_amr/rho_amr[0],'r+', label = "amrvac")
# plt.plot(y/H_eff, rho_non/rho_non[0], 'b-', label = "non-iso $\\rho$")
# plt.xlabel('$y/H_{eff}$', fontsize = 20)
# plt.ylabel('$\\rho/\\rho_0$', fontsize = 20)
# plt.legend()
#
# plt.figure(2)
# plt.title('Gas Energy $e_g$')
# plt.plot(y_amr,e_amr/e_amr[0],'r+', label = "amrvac")
# plt.plot(y/H_eff, e_non/e_non[0], 'b-', label = "non-iso $e_g$")
# plt.xlabel('$y/H_{eff}$', fontsize = 20)
# plt.ylabel('$e/e_0$', fontsize = 20)
# plt.legend()
#
# plt.figure(3)
# plt.title('Radiative Energy $E_r$')
# plt.plot(y_amr,r_e_amr/r_e_amr[0],'r+', label = "amrvac")
# plt.plot(y/H_eff, E_non/E_non[0], 'b-', label = "non-iso $E_r$")
# plt.xlabel('$y/H_{eff}$', fontsize = 20)
# plt.ylabel('$E_r/Er_0$', fontsize = 20)
# plt.legend()
#
# plt.figure(4)
# plt.title('Temperature $T$')
# plt.plot(tau_amr, T_amr/T_amr[0], 'r+', label = "amrvac")
# plt.plot(tau, T_non/T_non[0], 'b-', label = "non-iso $\\tau$")
# plt.xlabel('$\\tau$', fontsize = 20)
# plt.ylabel('$T/T_0$', fontsize = 20)
# plt.gca().invert_xaxis()
# plt.legend()

################################################################################

plt.figure(0)
plt.title('Gas Pressure')
plt.plot(y/H_eff, p_iso/p_bound, 'r--', label = "$p = p_0 e^{-y/H_{eff}}$")
plt.plot(y/H_eff, p_non/p_bound, 'b+', label = "$\\frac{dp_g}{dY} = \\frac{c p_g}{(ap_g + b){1/4}}$")
plt.xlabel('$y/H_{eff}$', fontsize = 20)
plt.ylabel('$p/p_0$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure(1)
plt.title('Density')
plt.plot(y/H_eff, rho_iso/rho_bound, 'r--', label = "$\\rho = p_g/a^2$")
plt.plot(y/H_eff, rho_non/rho_bound, 'b+', label = "$\\rho = \\frac{m_p \\mu}{k_b} \\frac{p_g}{(a p_g + b)^{1/4}}$")
plt.xlabel('$y/H_{eff}$', fontsize = 20)
plt.ylabel('$\\rho/\\rho_0$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure(2)
plt.title('Temperature')
plt.plot(y/H_eff,T_iso/T_bound, 'r--', label = "$T = \\frac{m_p \\mu}{k_b} a^2$")
plt.plot(y/H_eff,T_non/T_bound, 'b+', label = "$T = (a p_g + b)^{1/4}$")
plt.plot(y/H_eff,T_check/T_bound, 'bx', label = "$T = \\frac{3 F}{4 \\sigma}(\\tau - \\tau_0) + T_0$")
plt.plot(y/H_eff,T_check2/T_bound, 'b.', label = "$T = \\frac{m_p \\mu}{k_b} \\frac{p_g}{\\rho}$")
plt.xlabel('$y/H_{eff}$', fontsize = 20)
plt.ylabel('$T/T_0$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure(3)
plt.title('Temperature (Optical depth)')
plt.plot(tau,T_iso/T_bound, 'r--', label = "$T = \\frac{m_p \\mu}{k_b} a^2$")
plt.plot(tau,T_non/T_bound, 'b+', label = "$T = (a p_g + b)^{1/4}$")
plt.plot(tau,T_check/T_bound, 'bx', label = "$T = \\frac{3 F}{4 \\sigma}(\\tau - \\tau_0) + T_0$")
plt.plot(tau,T_check2/T_bound, 'b.', label = "$T = \\frac{m_p \\mu}{k_b} \\frac{p_g}{\\rho}$")
plt.gca().invert_xaxis()
plt.xlabel('$\\tau$', fontsize = 20)
plt.ylabel('$T/T_0$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure(4)
plt.title('Radiation energy')
plt.plot(y/H_eff, E_iso/E_iso[0], 'r--', label ='$E = \\frac{3 \\Gamma}{1-\\Gamma} p$')
plt.plot(y/H_eff, tau/tau_0,'k+', label = '$\\tau/\\tau_0$')
plt.plot(y/H_eff, E_non/E_non[0],'b+', label = '$E = \\frac{4 \\sigma}{c} T^4$')
plt.plot(y/H_eff, E_check/E_check[0],'bx', label = '$E = \\frac{4 \\sigma}{c} (a p_g + b)$')
plt.xlabel('$y/H_{eff}$', fontsize = 20)
plt.ylabel('$E/E_0, \\tau/ \\tau_0$', fontsize = 20)
plt.legend(loc = 'best')

plt.figure(5)
plt.title('Flux')
plt.plot(y[1:-1],F_check[1:-1]/F_check[1],'b+', label = '$c \lambda \\frac{\partial E}{\partial \\tau}$')
plt.plot(y[1:-1],F_non[1:-1]/F_non[1],'bx', label = '$\\frac{-c \lambda}{\kappa \\rho} \\nabla E$')
plt.plot(y[1:-1],np.ones(len(y))[1:-1],'r--', label = 'Flux')
plt.xlabel('$y/H_{eff}$', fontsize = 20)
plt.ylabel('$F/F_0$', fontsize = 20)
plt.legend(loc = 'best')

plt.show()
