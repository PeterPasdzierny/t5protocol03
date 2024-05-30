from matplotlib import pyplot as plt
import numpy as np

kb = 1.38e-23
ev_to_j = 1.6022e-19
j_to_ev = 6.242e18
kelvin = 273.15


def n_mol(T, emol):
    a = (2 * emol**0.5) / (np.pi**0.5 * (kb * T) ** 1.5)
    b = np.exp((-emol) / (kb * T))
    n_mol = a * b
    return n_mol


e_mol = np.linspace(0, 0.2, 101) * ev_to_j
temps = [-50 + kelvin, 0 + kelvin, 100 + kelvin]

plt.figure()
for t in temps:
    nmol = n_mol(t, e_mol)
    plt.plot(e_mol * j_to_ev, nmol, label=f"T = {t:.2f} K")
plt.xlim(0, 0.2)
plt.ylim(0, 1.6e20)
plt.grid()
plt.legend()

plt.show()
