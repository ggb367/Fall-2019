import math as m
import numpy as np
import numpy.linalg as lg
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import ode
import space_functions as sf


def JTwo(t, Y, mu, j, R):
    dY = np.zeros([6, 1])
    dY[0] = Y[3]
    dY[1] = Y[4]
    dY[2] = Y[5]
    r = np.sqrt(Y[0]**2+Y[1]**2+Y[2]**2)
    dY[3] = -(mu*Y[0]/r**3)*(1-j*1.5*(R/r)**2*(5*(Y[2]/r)**2-1))
    dY[4] = -(mu*Y[1]/r**3)*(1-j*1.5*(R/r)**2*(5*(Y[2]/r)**2-1))
    dY[5] = -(mu*Y[2]/r**3)*(1-j*1.5*(R/r)**2*(5*(Y[2]/r)**2-3))
    return dY

MU = 398600.4415
RE = 6378.1363
J2 = 0.00108248

def derivFcn(t, y):
    return JTwo(t, y, MU, J2, RE)

r_0 = [-5282.628, -4217.988, 296.511]
v_0 = [-4.541972, 5.885228, 2.043106]

E = sf.cart2elm(r_0, v_0, MU)
a = E[0]
P = 2*m.pi*np.sqrt(a**3/MU)
T0 = 0.0;
times = np.arange(0, 86400, 100)
tF = 86400
dT = 100
Y_0 = np.concatenate([r_0, v_0], axis=0)
rv = ode(derivFcn)

#  The integrator type 'dopri5' is the same as MATLAB's ode45()!
#  rtol and atol are the relative and absolute tolerances, respectively
rv.set_integrator('dopri5', rtol=1e-6, atol=1e-6)
rv.set_initial_value(Y_0, T0)
output = []
output.append(np.insert(Y_0, 0, T0))

# Run the integrator and populate output array with positions and velocities
while  rv.successful() and rv.t < tF:  # rv.successful() and
    rv.integrate(rv.t + dT)
    output.append(np.insert(rv.y, 0, rv.t))

#  Convert the output a numpy array for later use
output = np.array(output)
r = np.sqrt(np.power(output[:, 1], 2)+np.power(output[:, 2], 2)+np.power(output[:, 3], 2))
v = np.sqrt(np.power(output[:, 4], 2)+np.power(output[:, 5], 2)+np.power(output[:, 6], 2))
a = np.empty(np.size(r))
print("OUTPUTUTUTUTUTUT")
print(np.shape(output))
for i in range(np.size(r)):
    a[i] = np.divide(-MU, r[i]**2)
t = output[:, 0]
r_vec = np.empty([865, 3])
v_vec = np.empty([865, 3])
elms = np.empty([np.size(times)+1, 6])
for i in range(np.size(r)):
    r_vec[i, 0] = output[i, 1]
    r_vec[i, 1] = output[i, 2]
    r_vec[i, 2] = output[i, 3]
    v_vec[i, 0] = output[i, 4]
    v_vec[i, 1] = output[i, 5]
    v_vec[i, 2] = output[i, 6]
    elms[i] = sf.cart2elm(r_vec[i, :], v_vec[i, :], MU)
eng = .5*np.power(v, 2) - np.divide(MU, r)
plt.figure()

plt.subplot(3, 1, 1)
plt.plot(t/86400, r)
plt.ylabel('km', size=16)

plt.subplot(3, 1, 2)
plt.plot(t/86400, v)
plt.ylabel('km/s', size=16)

plt.subplot(3, 1, 3)
plt.plot(t/86400, a)
plt.ylabel('km/s^2', size=16)
plt.xlabel('Days', size=16)
plt.figure()
plt.plot(t/86400, eng)
plt.ylabel('Energy', size=16)
plt.xlabel('Days', size=16)

plt.figure()
plt.subplot(6, 1, 1)
plt.plot(t/3600, elms[:, 0])
plt.xlabel('Days')
plt.ylabel('a')
plt.subplot(6, 1, 2)
plt.plot(t / 3600, elms[:, 1])
plt.xlabel('Days')
plt.ylabel('e')
plt.subplot(6, 1, 3)
plt.plot(t / 3600, elms[:, 2])
plt.xlabel('Days')
plt.ylabel('i')
plt.subplot(6, 1, 4)
plt.plot(t / 3600, elms[:, 3])
plt.xlabel('Days')
plt.ylabel('RAAN')
plt.subplot(6, 1, 5)
plt.plot(t / 3600, elms[:, 4])
plt.xlabel('Days')
plt.ylabel('Arg of Periapsis')
plt.subplot(6, 1, 6)
plt.plot(t / 3600, elms[:, 5])
plt.xlabel('Hours')
plt.ylabel('Anomaly')
plt.show()