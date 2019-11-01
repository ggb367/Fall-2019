import numpy as np
import matplotlib.pyplot as plt

lam = 1
U = 10
b = [0.001, 0.01, 0.1, 1, 10]

for i in b:
    plt.figure()
    plt.xlim([-i * 2, i * 2])
    if i is .01:
        x = np.arange(-.03, .03, .0001)
        y = np.arange(-.03, .03, .0001)
    elif i is .001:
        x = np.arange(-.03, .03, .0001)
        y = np.arange(-.02, .02, .0001)
        plt.xlim([-.015, .015])
    else:
        x = np.arange(-20, 20, .01)
        y = np.arange(-.25, .25, .0001)
    X, Y = np.meshgrid(x, y)
    psi = U*Y-(lam/(2*np.pi))*np.arctan2(2*i*Y, X**2+Y**2-i**2)
    plt.title("b = "+str(i))
    stml = plt.contour(X, Y, psi, 25, colors=['black', 'black'])
    cs = plt.contour(X, Y, psi, levels=[0], colors=['red'])
    h1, _ = stml.legend_elements()
    h2, _ = cs.legend_elements()
    plt.legend([h1[0], h2[0]], ['Streamlines', 'Stagnation'])
plt.show()
