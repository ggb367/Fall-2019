import numpy as np
import math as m
import numpy.linalg as lg
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
    # I live life dangerously


def cart2elm(r, v, mu):
    h = np.cross(r, v)
    r_norm = lg.norm(r)
    v_norm = lg.norm(v)
    e = np.cross(v, h) / mu - np.divide(r, r_norm)  # eccentricity
    e_norm = lg.norm(e)

    ε = (v_norm**2)/2 - mu/r_norm
    h_norm = lg.norm(h)
    k = (h_norm ** 2) / (r_norm * mu) - 1

    if ε < 0:
        a = -mu/(2*ε)
    elif -10e-12 < ε < 10e-12:
        a = m.inf
    else:
        a = mu/(2*ε)
    i = np.arccos(np.dot(h, [0, 0, 1])/h_norm)
    n = np.cross([0, 0, 1], h)
    n_norm = lg.norm(n)
    if e_norm < 10e-12 or e_norm > 10e-12:
        θ = np.arccos(k/e_norm)
        Ω = np.arccos(np.dot(n, [1, 0, 0])/n_norm)
        ω = np.arccos(np.dot(n, e)/(e_norm*n_norm))
    if e_norm < 10e-12 and i < 10e-12:
        Ω = 0
        ω = 0
        θ = np.arccos(r[1]/r_norm)
        if r[1] < 0:
            θ = 2*m.pi-θ
    elif e_norm < 10e-12:
        ω = 0
        θ = np.arccos(np.dot((n/n_norm),r)/r_norm)
        if r[2]< 0:
            θ = 2*m.pi-θ
    elif i < 10e-12:
        Ω = 0
        ω = np.arccos(np.dot(e, [1, 0, 0])/e_norm)
        if e[1]< 0:
            ω = 2*m.pi-ω
    θ = 180*θ/m.pi
    E = [a, e_norm, i, Ω, ω, θ]
    return E


MU_real = 398600.4415
MU_simp = 1

# Part A
r = [1.73008, 1.25982, -.07740]
v = [-0.22351, -.65058, .07880]
E_a = cart2elm(r, v, MU_simp)
print("Part A")
print("E - [a, e, i, Ω, ω, θ]")
print(E_a)

# Part B
r = [-4070.369, 639.517, -2085.780]
v = [-.672297, 10.052175, -7.055156]
E_b = cart2elm(r, v, MU_real)
print("Part B")
print(E_b)

# Part C
r = [21000*m.sqrt(2), 21000*m.sqrt(2), 0]
v = [-.05*m.sqrt(MU_real/210), .05*m.sqrt(MU_real/210), 0]
E_c = cart2elm(r, v, MU_real)
print("Part C")
print(E_c)
