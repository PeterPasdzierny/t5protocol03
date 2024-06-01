from matplotlib import pyplot as plt
import numpy as np
from scipy import constants, integrate
import sympy as sp

kb = constants.Boltzmann
ev_to_j = constants.electron_volt
j_to_ev = 6.242e18
kelvin = 273.15
h = constants.Planck
c = constants.c
h_quer = constants.hbar
r = 116e-12
n_rots = [20,30]


def n_mol(T, emol):
    a = (2 * emol**0.5) / (np.pi**0.5 * (kb * T) ** 1.5)
    b = np.exp((-emol) / (kb * T))
    n_mol = a * b
    return n_mol

def n_mol_tilde(T, emol):
    a = np.exp((-emol) / (kb * T))
    f = sp.exp((-x) / (kb * T))
    integral = sp.integrate(f, (x, 0, sp.oo))#
    n_mol_tilde = a/integral
    return n_mol_tilde

def E_rot(n_rot):
    I = masses[0]*masses[1]/(masses[0]+masses[1])*r
    E =  (n_rot-h_quer)**2/(2*I)
    return E


x = sp.symbols('x')
f = sp.exp((-x) / (kb * 273.15))
integral = sp.integrate(f, (x, 0, sp.oo))


e_mol = np.linspace(0, 0.2, 101) * ev_to_j
temps = [-50 + kelvin, 0 + kelvin, 100 + kelvin]
masses = np.array([1.604e-2, 4.8e-2])/constants.Avogadro
speed = (2*e_mol/(masses[0]))**0.5




plt.figure()

for t in temps:
    nmol_T = n_mol(t, e_mol)
    plt.plot(e_mol * j_to_ev, nmol_T, label=f"T = {t-kelvin:.0f} °C") 
    

plt.xlim(0, 0.2)
plt.ylim(0, 1.6e20)
plt.xlabel('$E_{trl}$ [eV]')
plt.ylabel('$N_{mol}$')
plt.grid()
plt.legend()




linestyles = ['solid', 'dashed']
colors = ['tab:blue', 'tab:orange', 'tab:green']
mass_names = ['Methane', 'Ozone']



plt.figure()

for i,m in enumerate(masses):
    e_mol =  m/2 * speed**2 
    
    for j, t in enumerate(temps):
        if i == 0:
            label = f"T = {t-kelvin:.0f} °C"
        else:
            label = None
        nmol_v = n_mol(t, e_mol)
        
        plt.plot(speed, nmol_v, label=label, linestyle = linestyles[i], c= colors[j]) 
    
        print(f'Most probable velocity for {mass_names[i]} is {speed[np.argmax(nmol_v)]} ms⁻¹')

        e_diff = (m/2 * 615**2)- e_mol[[np.argmax(nmol_v)]]
        lam = (c)*h/e_diff


plt.xlabel('$v$ [ms⁻¹]')
plt.ylabel('$N_{mol}$')
plt.grid()
plt.legend()



print(f'The ozon photon will have an energy of {lam*1e9} nm.')



plt.figure()

for t in temps:
    nmol_T = n_mol_tilde(t, e_mol)
    plt.plot(e_mol * j_to_ev, nmol_T, label=f"T = {t-kelvin:.0f} °C") 

plt.xlabel('$E$ [eV]')
plt.ylabel('$N_{mol}$')
plt.xscale('log')
plt.grid()
plt.legend()
#plt.show()




for i, n_rot in enumerate(n_rots):
    
    E = E_rot(n_rot)
    
    nmol_T_dif = n_mol_tilde(273.15, E)- n_mol_tilde(373.15, E)
   
    print(f'The difference for {mass_names[i]} is {nmol_T_dif}.')
        

