import matplotlib.pyplot as plt
import numpy as np

du = np.linspace(1, 5, 100)
g = np.divide(1, np.square(du))
plt.plot(du, g)
plt.xlabel('Distance in Earth Radii')
plt.ylabel('g/g_0')
plt.title('Inverse Square Decay of Gravity Force')
plt.show()
