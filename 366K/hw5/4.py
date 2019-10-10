import numpy as np
import math as m


def elm2cart(E, mu):
    # E - [a, e, i, Ω, ω, θ]
    a = E[0]
    e = E[1]
    i = m.pi * E[2] / 180
    Ω = m.pi * E[3] / 180
    ω = m.pi * E[4] / 180
    θ = m.pi * E[5] / 180
    p = a*(1 - e**2)
    r_pqw = np.array([(p/(1+e*np.cos(θ)))*np.cos(θ), (p/(1+e*np.cos(θ)))*np.sin(θ), 0])
    v_pqw = np.array([np.sqrt(mu/p)*(-np.sin(θ)), np.sqrt(mu/p)*(e+np.cos(θ)), 0])
    # R_3(-Ω)R_1(-i)R_3(-ω)
    c1 = np.cos(-ω)
    c2 = np.cos(-i)
    c3 = np.cos(-Ω)
    s1 = np.sin(-ω)
    s2 = np.sin(-i)
    s3 = np.sin(-Ω)
    q1 = np.array([c1*c3-c2*s1*s3, c3*s1+c1*c2*s3, s3*s2])
    q2 = np.array([-c1*s3-c3*c2*s1, c1*c2*c3-s1*s3, c1*s2])
    q3 = np.array([s1*s2, -c1*s2, c2])
    Q = np.array([q1, q2, q3])
    r = np.matmul(Q, r_pqw)
    v = np.matmul(Q, v_pqw)
    return r, v


MU_real = 398600.4415
MU_simp = 1
 #2
E = [1, .82, 171, 22, 130, 216]
r_2, v_2 = elm2cart(E, MU_simp)
print("Part 2:")
print("r = " + str(r_2))
print("v = " + str(v_2))
# Part A
E = [26778.9, 0.40000, 142.43, 177.92, 105.30, 162.17]
r_a, v_a = elm2cart(E, MU_real)
print("Part A:")
print("r = " + str(r_a))
print("v = " + str(v_a))

# Part B
E = [1.4, 0.55069, 117.14, 268.84, 177.95, 188.18]
r_b, v_b = elm2cart(E, MU_simp)
print("Part B:")
print("r = " + str(r_b))
print("v = " + str(v_b))
