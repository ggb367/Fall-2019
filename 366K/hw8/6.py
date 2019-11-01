import numpy as np
import matplotlib.pyplot as plt

MU = 1
RE = 1
JDos = -.01

i = np.linspace(0, np.pi, 180)
e = [0.1, 0.3, 0.5]
a = np.empty([np.size(i), np.size(e)])
for k in range(np.size(e)):
    a[:, k] = np.multiply(-(1.5*(np.sqrt(MU)*JDos*RE**2)/((1-e[k]**2)**2*3.162e-8)), np.cos(i))**(2/7)
    for l in a[:, k]:
        index = np.where(a[:, k] == l)
        rp = l*(1-e[k])
        if rp < RE:
            a[index, k] = None
plt.plot(i*180/np.pi, a)
#plt.ylim([RE, 30])
plt.xlabel("Inclination [deg]")
plt.ylabel("Semimajor Axis [DU]")
plt.title("Semimajor Axis vs. Inclination")
plt.show()
