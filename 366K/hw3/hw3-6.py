import numpy as np
import numpy.linalg as lg

MU = 1  # DU^3/TU^2

r = [-1.2, 1.5, 0]
v = [0, 0, .922]
h = np.cross(r, v)
print("The momentum vector is: "+str(h) + "DU^2/TU")
v_norm = lg.norm(v)
r_norm = lg.norm(r)
eps = (v_norm*v_norm)/2 - MU/r_norm
print("Total Energy is: "+str(eps))
e = np.cross(v, h)/MU - np.divide(r, r_norm)
print("Eccentricity vector is: " + str(e))
e_norm = lg.norm(e)
h_norm = lg.norm(h)
k = (h_norm*h_norm)/r_norm - 1
theta = np.arccos(k/e_norm)
if np.dot(r, v) < 0:
    theta = 2*np.pi() - theta
print("Anomaly is: " + str(theta))
print("The magnitude of the eccentricity is: " + str(e_norm))
a = (r_norm*(1+e_norm*np.cos(theta)))/(1-e_norm*e_norm)
print("Semi-major axis length: " + str(a))
r_p = a*(1-e_norm)
print("Periapsis length: " + str(r_p))
r_a = a*(1+e_norm)
print("Apoapsis length: " + str(r_a))
