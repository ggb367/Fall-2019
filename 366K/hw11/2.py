import matplotlib.pyplot as plt
import numpy as np

import space_functions as sf

mu = 0.0121506037932213

tf = 15
dT = 0.001
lagrange = sf.lagrange(mu)  # calculate lagrange points from mu
r_0 = [-.08, -.03, 0]  #this will change depending on the initial conds
v_0 = [3.5, -3.2, 0]
r_vec, _ = sf.CRTBP_prop_rk(r_0, v_0, 0, tf, dT, mu)  # propogate in CRTBP Frame
r_vec = np.multiply(384.400, r_vec)  # transform the position vector to SI Units
time_series = np.arange(0, 15, .001)
plt.figure()  # Synodic Frame
for i in [0, 1, 2, 3, 4]:  # I am ashamed of this for loop, but what is done is done
    if i is 0:  # plot lagrange points
        plt.scatter(lagrange[i, 0]*384.400, lagrange[i, 1]*384.400, color='r', label='Lagrange Points')
    plt.scatter(lagrange[i, 0] * 384.400, lagrange[i, 1] * 384.400, color='r')
plt.plot(r_vec[:, 0], r_vec[:, 1], label='Trajectory')  # plot trajectory
plt.scatter(384.400, 0, label='Secondary Body')
plt.scatter(0, 0, label='Primary Body')
plt.title("Orbit in Synodic Frame")
plt.legend()
plt.xlabel("X [x1000 km]")
plt.ylabel("Y [x1000 km]")

plt.figure()  # Inertial Frame
r_inrt = np.empty(np.shape(r_vec))
r_secondary = np.empty([np.size(time_series), 2])
for i in range(np.size(r_vec, 0)-2):
    Q_synodic = sf.R3(-time_series[i])  # calculate transform
    r_inrt[i, :] = np.matmul(Q_synodic, r_vec[i, :]-np.multiply(-mu,[1, 0, 0]))  # transform to inertial frame
    r_secondary[i, :] = [384.400*np.cos(time_series[i]), 384.400*np.sin(time_series[i])]  # calculate secondary orbit
plt.scatter(0, 0, label="Primary Body")  # TODO Fix zero value artifact
plt.plot(r_inrt[:, 0], r_inrt[:, 1], label="Trajectory")
plt.plot(r_secondary[:,0], r_secondary[:, 1], label="Secondary Trajectory")
plt.legend()
plt.title("Orbit in Earth Centered Inertial Frame")
plt.xlabel("X [x1000 km]")
plt.ylabel("Y [x1000 km]")
plt.show()


