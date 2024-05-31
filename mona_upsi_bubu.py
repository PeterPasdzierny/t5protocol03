from matplotlib import pyplot as plt
import numpy as np
from scipy import constants 

kb = constants.Boltzmann
ev_to_j = constants.electron_volt
j_to_ev = 6.242e18
kelvin = 273.15
h = constants.Planck
c = constants.c



def n_mol(T, emol):
    a = (2 * emol**0.5) / (np.pi**0.5 * (kb * T) ** 1.5)
    b = np.exp((-emol) / (kb * T))
    n_mol = a * b
    return n_mol


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
plt.show()



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

plt.show()

print(f'The ozon photon will have an energy of {lam*1e9} nm.')
