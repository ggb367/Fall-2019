import numpy as np
import matplotlib.pyplot as plt

lam = 1
U = 10
source = .1
sink = 1
b = 1


x = np.linspace(-1.5, 1.5, 100)
y = np.linspace(-0.5, 0.5, 100)
X, Y = np.meshgrid(x, y)
psi = U * Y + source/(2*np.pi)*np.arctan2(Y, X+b)-sink/(2*np.pi)*np.arctan2(Y, X+b)
stml = plt.contour(X, Y, psi, 30)
plt.title("b = 1, Sink Stronger")
cs = plt.contour(X, Y, psi, levels=[0], colors=['red'])
h1, _ = stml.legend_elements()
h2, _ = cs.legend_elements()
plt.legend([h1[0], h2[0]], ['Streamlines', 'Stagnation'])

plt.show()
