import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt

MU = 1  # DU^3/TU^2

# initial conditions
r = [-1.2, 1.5, 0]
v = [0, 0, .922]
theta = np.linspace(0, 2*np.pi, 361)  # numpy uses radians
# calculate what we need to find the velocity and radius
h = np.cross(r, v)  # momentum
e = np.cross(v, h)/MU - np.divide(r, lg.norm(r))  # eccentricity
e_norm = lg.norm(e)
r = lg.norm(r)
v = lg.norm(v)
cos = np.cos(theta)
h = lg.norm(h)
# calculate r and v as a function of the anomaly
r = np.multiply((h*h/MU), np.divide(1, (1+np.multiply(e_norm, cos))))  # r as a func of theta
a = (r[205]*(1+e_norm*np.cos(theta[205])))/(1-e_norm*e_norm)  # use any r and theta paring to find a
v = np.sqrt(MU*(np.divide(2, r))-1/a)  # use r and a to find v
# Find the energies
KE = 0.5*np.square(v)
PE = np.divide(MU, r)
TE = KE-PE
dTE = TE - (KE[0]-PE[0])
# plot everything
plt.plot(theta, KE, label="Kinetic Energy")
plt.plot(theta, PE, label="Potential Energy")
plt.plot(theta, dTE, label="Change in Total Energy")
plt.legend()
plt.xlabel('Anomaly(θ)')
plt.ylabel('Energy')
plt.title('Energy vs Anomaly')
plt.show()
