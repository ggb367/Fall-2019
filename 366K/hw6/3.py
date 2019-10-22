import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as lg
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
jdos = 0.00108248

def derivFcn(t, y):
    return JTwo(t, y, MU, jdos, RE)

r_0 = [-5282.628, -4217.988, 296.511]
v_0 = [-4.541972, 5.885228, 2.043106]

T0 = 0.0;
times = np.arange(0, 86400, 100)
tF = 86400
dT = 100
Y_0 = np.concatenate([r_0, v_0], axis=0)
rv = ode(derivFcn)

#  The integrator type 'dopri5' is the same as MATLAB's ode45()!
#  rtol and atol are the relative and absolute tolerances, respectively
rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
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
r_norm = lg.norm(r)
t = output[:, 0]

r_vec = np.empty([865, 3])
v_vec = np.empty([865, 3])
a = np.empty(np.size(r))
U = np.empty(np.size(r))
eng = np.empty(np.size(r))
elms = np.empty([np.size(times)+1, 6])
eng_0 = (np.dot(v_0, v_0))/2 - MU/r_norm-(MU/r_norm)*(jdos/2)*(RE/r_norm)**2*(3*(r[2]/r_norm)**2-1)

for i in range(np.size(r)):
    r_vec[i, 0] = output[i, 1]
    r_vec[i, 1] = output[i, 2]
    r_vec[i, 2] = output[i, 3]
    v_vec[i, 0] = output[i, 4]
    v_vec[i, 1] = output[i, 5]
    v_vec[i, 2] = output[i, 6]
    Y = [float(output[i, 1]), float(output[i, 2]), float(output[i, 3]), float(output[i, 4]), float(output[i, 5]), float(output[i, 6])]
    dY = JTwo(t[i], Y, MU, jdos, RE)
    a[i] = -lg.norm([dY[3], dY[4], dY[5]])
    elms[i] = sf.cart2elm(r_vec[i, :], v_vec[i, :], MU)
    r_vec_norm = lg.norm(r_vec[i, :])
    U[i] = MU/r_vec_norm-(MU/r_vec_norm)*(jdos/2)*(RE/r_vec_norm)**2*(3*(r_vec[i, 2]/r_vec_norm)**2-1)
    eng[i] = (np.dot(v_vec[i, :], v_vec[i, :])/2 - U[i]) - eng_0
E_0 = sf.cart2elm(r_0, v_0, MU, deg=False)
OMG_dot = (-(1.5*(np.sqrt(MU)*jdos*RE**2)/(((1-E_0[1]**2)**2)*E_0[0]**3.5)))*np.cos(E_0[2])
print("Omega dot: "+ str(OMG_dot)+" rad/s")

plt.figure()
plt.suptitle("Position, Velocity, & Acceleration vs. Time")
plt.subplot(3, 1, 1)
plt.plot(t/3600, r)
plt.ylabel('km', size=16)
plt.subplot(3, 1, 2)
plt.plot(t/3600, v)
plt.ylabel('km/s', size=16)
plt.subplot(3, 1, 3)
plt.plot(t/3600, a)
plt.ylabel('$km/s^2$', size=16)
plt.xlabel('Hours', size=16)

plt.figure()
plt.title("Change in Energy vs. Time")
plt.plot(t/3600, eng)
plt.ylabel('Energy', size=16)
plt.xlabel('Hours', size=16)

plt.figure()
plt.suptitle("Orbital Elements vs. Time")
plt.subplot(3, 2, 1)
plt.plot(t/3600, elms[:, 0])
plt.ylabel('a [km]')
plt.subplot(3, 2, 2)
plt.plot(t / 3600, elms[:, 1])
plt.ylabel('e')
plt.subplot(3, 2, 3)
plt.plot(t / 3600, elms[:, 2])
plt.ylabel('i [deg]')
plt.subplot(3, 2, 4)
plt.plot(t / 3600, elms[:, 3])
plt.ylabel('RAAN [deg]')
plt.subplot(3, 2, 5)
plt.plot(t / 3600, elms[:, 4])
plt.ylabel('Arg of Periapsis [deg]')
plt.subplot(3, 2, 6)
plt.plot(t / 3600, elms[:, 5])
plt.xlabel('Hours')
plt.ylabel('Anomaly [deg]')
plt.show()