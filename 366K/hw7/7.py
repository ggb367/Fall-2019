import numpy as np
import matplotlib.pyplot as plt

MU = 398600.4415
RE = 6378.1363
JDos = 0.00108248
f = 3.353e-3

i = np.linspace(95*np.pi/180, np.pi, 180-95)
e = [0.1, 0.3, 0.4, 0.5, 0.9]
a = np.empty([np.size(i), np.size(e)])
for k in range(np.size(e)):
    a[:, k] = np.multiply(-(1.5*(np.sqrt(MU)*JDos*RE**2)/((1-e[k]**2)**2*1.991e-7)), np.cos(i))**(2/7)
plt.plot(i*180/np.pi, a)
plt.xlabel("Inclination [deg]")
plt.ylabel("Semimajor Axis [km]")
plt.title("Semimajor Axis vs. Inclination")
plt.show()
