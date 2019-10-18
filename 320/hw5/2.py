import matplotlib.pyplot as plt
import numpy as np

c = [1, 2, 3, 4, 5]
circle = []
fig, ax = plt.subplots()
for i in range(np.size(c)):
    circle.append(plt.Circle((0, 0), c[i], fill=False))
    ax.add_artist(circle[i])
ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
plt.show()
