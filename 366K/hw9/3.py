import space_functions as sf
import numpy as np
import numpy.linalg as lg
import matplotlib.pyplot as plt
import scipy.io as spio

MU = 398600
omega_earth = 2*np.pi/86164
f = 3.35e-3

# Orbit A
init_cond = [7000, 0.01, 45, 0, -65, 0]
theta_g_i = 0
t_0 = 0
P = 2*np.pi*np.sqrt(init_cond[0]**3/MU)  # seconds
[r_a, v_a] = sf.elm2cart(init_cond, MU)
[r, v] = sf.orbit_prop_rk(r_a, v_a, t_0, P*3, 60)
time_series = np.arange(0, 3*P, 60)
step = 0
phi_gd = np.empty(np.shape(time_series))
lam = np.empty(np.shape(time_series))
for t_n in time_series:
    theta_g = theta_g_i+omega_earth*(t_n-t_0)
    r_ecf = np.matmul(sf.R3(theta_g), r[step, :])
    r_norm = lg.norm(r_ecf)
    phi_gc = np.arcsin(r_ecf[2]/r_norm)
    lam[step] = np.arctan2(r_ecf[1], r_ecf[0])
    phi_gd[step] = np.arctan2(np.tan(phi_gc), (1-f)**2)
    step = step+1
# Convert to degrees
lam = np.multiply(lam, 180/np.pi)
phi_gd = np.multiply(phi_gd, 180/np.pi)
#  Set the name of the file with the coastline data
coastFile = 'earth_coastline.mat'
#  Load the file.
earth_coastline = spio.loadmat(coastFile)['earth_coastline']
#  Generate the figure for the plot
plt.figure()
#  Plot the coastlines
plt.plot(earth_coastline[:, 0], earth_coastline[:, 1], 'k')
#  Set the correct aspect ratio for the plot
plt.axis('square')
plt.axis([-180, 180, -90, 90])
plt.scatter(lam, phi_gd, marker='.', c='r')
#  Show the plot
plt.show()
# Orbit B
