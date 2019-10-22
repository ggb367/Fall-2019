import numpy as np
import matplotlib.pyplot as plt

MU = 398600.4415
RE = 6378.1363
JDos = 0.00108248
f = 3.353e-3

i = 116.565
a = np.linspace(RE, 35000, 100000)
e = np.sqrt(1+np.multiply(1.5*(np.sqrt(MU)*JDos*RE**2)/1.991e-7*np.cos(i*np.pi/180), a**(-7/2)))
plt.plot(e, a)
plt.xlabel("Eccentricity")
plt.ylabel("Semimajor Axis [km]")
plt.title("a vs. e for a Critically Inclined, Sun-Synchronous Orbit")
plt.show()
