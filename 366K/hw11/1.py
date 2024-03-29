import numpy as np
from scipy.optimize import fsolve

# givens
m1 = 6e24;
m2 = 3e24;
#normalize
mu = m2/(m1+m2)  # calculate mu
f = lambda r_x : r_x - (1-mu)*(r_x+mu)/np.abs(r_x+mu)**3 - mu*(r_x-(1-mu))/np.abs(r_x+mu-1)**3  # function for y = 0

r_x = np.array([-2, 0, 2])
r_0_roots = np.array(fsolve(f, r_x))  # calculate y=0 points
roots_x = np.append(r_0_roots, [.5-mu, .5-mu])  # calculate y!=0 points
roots_y = np.array([0, 0, 0, np.sqrt(3)/2, -np.sqrt(3)/2])
roots = np.column_stack((roots_x, roots_y))
print(roots)
