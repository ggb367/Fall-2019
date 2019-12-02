import numpy as np
import space_functions as sf
from scipy.optimize import fsolve

m1 = 6e24;
m2 = 3e24;
#normalize
mu = m2/(m1+m2)
f = lambda r_x : r_x - (1-mu)*(r_x+mu)/np.abs(r_x+mu)**3 - mu*(r_x-(1-mu))/np.abs(r_x+mu-1)**3

r_x = np.array([-2, 0, 2])
r_0_roots = np.array(fsolve(f, r_x))
roots_x = np.append(r_0_roots, [.5-mu, .5-mu])
roots_y = np.array([0, 0, 0, np.sqrt(3)/2, -np.sqrt(3)/2])
roots = np.column_stack((roots_x, roots_y))
print(roots)
