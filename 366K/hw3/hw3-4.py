import numpy as np
import matplotlib.pyplot as plt

theta = np.linspace(0, 2*np.pi, 361)  # all the radians
e = np.linspace(0, .8, 5)
for i in range(np.size(e)):
    e_temp = e[i]
    r = 10000*np.divide(1, (1+np.multiply(e_temp, np.cos(theta))))
    x = np.multiply(r, np.cos(theta))
    y = np.multiply(r, np.sin(theta))
    plt.plot(x, y)
plt.show()
