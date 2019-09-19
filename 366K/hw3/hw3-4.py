import numpy as np
import matplotlib.pyplot as plt

R_P = 10000  # periapsis
R_E = 6300  # radius of the earth

theta = np.linspace(0, 2*np.pi, 361)  # all the radians
e = [0, .2, .4, .6, .8]
for i in range(np.size(e)):
    e_temp = e[i]
    # find radius
    r = (1+e_temp)*R_P*np.divide(1, (1+np.multiply(e_temp, np.cos(theta))))
    r = np.divide(r, R_E)
    # convert from polar to cartesian
    x = np.multiply(r, np.cos(theta))  # x = rcos(theta)
    y = np.multiply(r, np.sin(theta))  # y = rsin(theta)
    plt.plot(x, y, label=str(e_temp))
plt.grid()
#plt.plot(0, 0, 'bo')
circ = plt.Circle((0, 0), 1, color='b')
plt.gcf().gca().add_artist(circ)
plt.legend(bbox_to_anchor=(1.165, 1), loc="upper right", borderaxespad=0.)
plt.title("Orbits with Varying Eccentricity")
plt.axis('equal')
plt.xlabel('Distance [R⊕]')
plt.ylabel('Distance [R⊕]')
plt.show()
