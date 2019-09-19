import numpy as np
import matplotlib.pyplot as plt

R_P = 10000

theta = np.linspace(0, 2*np.pi, 361)  # all the radians
e = [0, .2, .4, .6, .8]
for i in range(np.size(e)):
    e_temp = e[i]
    # find radius
    r = R_P*np.divide(1, (1+np.multiply(e_temp, np.cos(theta))))
    r = np.divide(r, R_P)
    # convert from polar to cartesian
    x = np.multiply(r, np.cos(theta))  # x = rcos(theta)
    y = np.multiply(r, np.sin(theta))  # y = rsin(theta)
    plt.plot(x, y, label=str(e_temp))
plt.grid()
plt.plot(0, 0, 'bo')
plt.legend(bbox_to_anchor=(1.2, 1), loc="upper right", borderaxespad=0.)
plt.title("Orbits with Varying Eccentricity")
plt.axis('equal')
plt.show()
