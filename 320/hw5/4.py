import matplotlib.pyplot as plt
import numpy as np

c = np.arange(13)
print(c)
x = np.linspace(-5, 5, 1000)
fig, ax = plt.subplots()
y_eq = np.empty([1000, np.size(c)*2])
y_v = np.empty([1000, np.size(c)*2])
x1, y = np.meshgrid(x, x)
for i in range(np.size(c)):
    y_v[:, -i+np.size(c)] = np.multiply(c[i], np.divide(1, x))
    y_v[:, i] = np.multiply(-c[i], np.divide(1, x))
    plt.contour(x1, y, (x ** 2 - y**2 - c[i]), [1])
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
plt.plot(x, y_v)
plt.show()
