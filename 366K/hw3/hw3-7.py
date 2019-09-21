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
cos = np.cos(theta)
h = lg.norm(h)
# calculate r and v as a function of the anomaly
r = np.multiply((h*h/MU), np.divide(1, (1+np.multiply(e_norm, cos))))  # orbit equation
a = (r[205]*(1+e_norm*np.cos(theta[205])))/(1-e_norm*e_norm)  # use any r and theta paring to find a
v = np.sqrt(MU*(np.divide(2, r))-1/a)  # use r and a to find v, from vis-viva
# Find the energies
KE = 0.5*np.square(v)
PE = -np.divide(MU, r)
TE = KE+PE
dTE = TE - (KE[0]+PE[0])
# plot everything
fig, axs = plt.subplots(3, sharex=True, gridspec_kw={'hspace': .3})
axs[0].plot(theta, KE)
axs[0].set_title("Kinetic Energy")
axs[1].plot(theta, PE)
axs[1].set_title("Potential Energy")
axs[2].plot(theta, dTE)
axs[2].set_title("Change in Total Energy")
axs[2].set_ylim([-0.5, 0.5])
plt.xlabel('Anomaly (θ)')
plt.ylabel('Energy [DU^2/TU^2]')
fig.suptitle('Energy vs Anomaly', fontsize=16, fontweight='bold')
plt.show()
