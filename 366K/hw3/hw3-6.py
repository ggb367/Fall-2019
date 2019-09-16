import numpy as np
import numpy.linalg as lg

MU = 1  # DU^3/TU^2

r = [-1.2, 1.5, 0]
v = [0, 0, .922]
h = np.cross(r, v)
print("H IS "+str(h))
v_norm = lg.norm(v)
r_norm = lg.norm(r)
eps = (v_norm*v_norm)/2 - MU/r_norm
print("EPS IS "+str(eps))
e = np.cross(v, h)/MU - np.divide(r, r_norm)
print("E IS " + str(e))
e_norm = lg.norm(e)
h_norm = lg.norm(h)
k = (h_norm*h_norm)/r_norm - 1
theta = np.arccos(k/e_norm)
if np.dot(r, v) < 0:
    theta = 2*np.pi() - theta
print("THETA IS " + str(theta))
print("E IS " + str(e_norm))
a = (r_norm*(1+e_norm*np.cos(theta)))/(1-e_norm*e_norm)
print("A IS " + str(a))
r_p = a*(1-e_norm)
print("R_P IS " + str(r_p))
r_a = a*(1+e_norm)
print("R_A IS " + str(r_a))
