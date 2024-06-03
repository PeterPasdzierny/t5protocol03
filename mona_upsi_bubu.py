from matplotlib import pyplot as plt
import numpy as np
from scipy import constants, integrate
import sympy as sp

kb = constants.Boltzmann
ev_to_j = constants.electron_volt
e = constants.elementary_charge
j_to_ev = 6.242e18
kelvin = 273.15
h = constants.Planck
n = constants.Avogadro
c = constants.c
h_quer = constants.hbar
r = 116e-12
n_rots = [20, 30]
me = constants.m_e
eps0 = 8.854e-12


def n_mol(T, emol):
    a = (2 * emol**0.5) / (np.pi**0.5 * (kb * T) ** 1.5)
    b = np.exp((-emol) / (kb * T))
    n_mol = a * b
    return n_mol


def n_mol_v(T, m, v):
    a = (m / (2 * np.pi * kb * T)) ** 1.5 * 4 * np.pi * v**2
    b = np.exp((-m * v**2) / (2 * kb * T))
    n_mol_v = a * b
    return n_mol_v


def get_v_max(T, m):
    return np.sqrt(2 * kb * T / m)


def n_mol_disc(T, emol):
    n_mol_disc = np.exp((-emol) / (kb * T)) / integrate.trapezoid(
        np.exp(-emol / (kb * T)), dx=0.1
    )
    return n_mol_disc


def E_orb(n_orb):
    return -me * e**4 / (8 * eps0**2 * h**2 * n_orb**2)


def E_rot(n_rot, I):
    return n_rot * (n_rot + 1) * h**2 / (8 * np.pi**2 * I)


def E_vib(n_vib, nu_norm):
    return (n_vib + 0.5) * h * nu_norm


# def n_mol_tilde(T, emol):
#     a = np.exp((-emol) / (kb * T))
#     f = sp.exp((-x) / (kb * T))
#     integral = sp.integrate(f, (x, 0, sp.oo))  #
#     n_mol_tilde = a / integral
#     return n_mol_tilde
#
#
# def E_rot(n_rot):
#     I = masses[0] * masses[1] / (masses[0] + masses[1]) * r
#     E = (n_rot - h_quer) ** 2 / (2 * I)
#     return E


# x = sp.symbols("x")
# f = sp.exp((-x) / (kb * 273.15))
# integral = sp.integrate(f, (x, 0, sp.oo))


# 1a
e_mol = np.linspace(0, 0.2, 1001) * ev_to_j
temps = [-50 + kelvin, 0 + kelvin, 100 + kelvin]

plt.figure()

for t in temps:
    nmol_T = n_mol(t, e_mol)
    plt.plot(e_mol * j_to_ev, nmol_T, label=f"T = {t-kelvin:.0f} °C")

plt.xlim(0, 0.2)
plt.ylim(0, 1.6e20)
plt.xlabel("$E_{trl}$ [eV]")
plt.ylabel("$N_{mol}$")
plt.grid()
plt.legend()

############################################################

# 1b)
vels = np.linspace(0, 1500, 10001)
masses = np.array([16.04e-3, 48e-3]) / n

# linestyles = ["solid", "dashed"]
# colors = ["tab:blue", "tab:orange", "tab:green"]
mass_names = ["Methane", "Ozone"]
#
plt.figure()

for t in temps:
    for m in masses:
        nmol_T = n_mol_v(t, m, vels)
        plt.plot(vels, nmol_T, label=f"T = {t-kelvin:.0f} °C")
        print(
            f"Most probable v for {list(mass_names)[list(masses).index(m)]} @ {t-kelvin} °C {get_v_max(t, m)}"
        )

plt.xlim(0, 1500)
# plt.ylim(0, 1.6e20)
plt.xlabel("$v$ [m/s]")
plt.ylabel("$N_{mol}$")
plt.grid()
plt.legend()

# 1c
e_diff = masses[1] / 2 * (615**2 - get_v_max(temps[2], masses[1]) ** 2)
lam = (c) * h / e_diff
print(f"The ozon photon will have a wavelength of {lam*1e6} um.")

#######################################################################

# 2a
n_rot = np.linspace(0, 500, 1001)
n_vib = np.linspace(0, 60, 1001)
n_orb = np.linspace(1, 8, 1001)
e_mol_rot = E_rot(n_rot, I=7.1501e-46)
e_mol_vib = E_vib(n_vib, nu_norm=136600 * c)
e_mol_orb = E_orb(n_orb)

plt.figure()

for t in temps:
    nmol_disc = n_mol_disc(t, e_mol_rot)
    plt.plot(e_mol_rot * j_to_ev, nmol_disc, label=f"T = {t-kelvin:.0f} °C", c="b")
    nmol_disc = n_mol_disc(t, e_mol_vib)
    plt.plot(e_mol_vib * j_to_ev, nmol_disc, label=f"T = {t-kelvin:.0f} °C", c="g")
    nmol_disc = n_mol_disc(t, e_mol_orb)
    plt.plot(
        e_mol_orb * j_to_ev + 13.6, nmol_disc, label=f"T = {t-kelvin:.0f} °C", c="r"
    )

plt.xlabel("$E$ [eV]")
plt.ylabel("$N_{mol}$")
plt.xscale("log")
plt.yscale("log")
plt.grid()
plt.legend()


# for i, n_rot in enumerate(n_rots):
#     E = E_rot(n_rot)
#
#     nmol_T_dif = n_mol_tilde(273.15, E) - n_mol_tilde(373.15, E)
#
#     print(f"The difference for {mass_names[i]} is {nmol_T_dif}.")


plt.tight_layout()
plt.show()
