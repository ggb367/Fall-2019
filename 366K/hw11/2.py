import space_functions as sf
import numpy as np
import matplotlib.pyplot as plt

mu = 0.0121506037932213

tf = 3
lagrange = sf.lagrange(mu)
r_0 = [lagrange[1, 0], lagrange[1, 1], 0]
v_0 = [0, 0, 0]
r_vec, v_vec = sf.CRTBP_prop_rk(r_0, v_0, 0, tf, 0.001, mu)
r_vec = np.multiply(384.400, r_vec)
time_series = np.arange(0, 15, .001)
plt.figure()
plt.plot(r_vec[:, 0], r_vec[:, 1], label='Trajectory')
plt.scatter(384.400, 0, label='Secondary Body')
plt.scatter(0, 0, label='Primary Body')
for i in [0, 1, 2, 3, 4]:
    if i is 0:
        plt.scatter(lagrange[i, 0]*384.400, lagrange[i, 1]*384.400, color='r', label='Lagrange Points')
    plt.scatter(lagrange[i, 0] * 384.400, lagrange[i, 1] * 384.400, color='r')
r_inrt = np.empty(np.shape(r_vec))
r_secondary = np.empty([np.size(time_series), 2])
for i in range(np.size(r_vec, 0)-2):
    Q_synodic = sf.R3(-time_series[i])
    r_inrt[i, :] = np.matmul(Q_synodic, r_vec[i, :]-np.multiply(-mu,[1, 0, 0]))
    r_secondary[i, :] = [384.400*np.cos(time_series[i]), 384.400*np.sin(time_series[i])]
plt.title("Orbit in Synodic Frame")
plt.legend()
plt.xlabel("X [x1000 km]")
plt.ylabel("Y [x1000 km]")
plt.figure()
plt.scatter(0, 0, label="Primary Body")
plt.plot(r_inrt[:, 0], r_inrt[:, 1], label="Trajectory")
plt.plot(r_secondary[:,0], r_secondary[:, 1], label="Secondary Trajectory")
plt.legend()
plt.title("Orbit in Earth Centered Inertial Frame")
plt.xlabel("X [x1000 km]")
plt.ylabel("Y [x1000 km]")
plt.show()


