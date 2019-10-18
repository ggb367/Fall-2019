import math as m
import numpy as np
import numpy.linalg as lg
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.integrate import ode
import space_functions as sf

MU = 398600.4415
RE = 6378.1363
J2 = 0.00108248
#r_0_norm = lg.norm(r_0)
#v_0_norm = lg.norm(v_0)
#v_r0 = np.dot(r_0, v_0)/r_0_norm
#v_r0_norm = lg.norm(v_r0)
#h = r_0_norm*np.sqrt(v_0_norm**2-v_r0_norm**2)


def two_body_orbit(t, Y, mu):
    dY = np.empty([6, 1])
    dY[0] = Y[3]
    dY[1] = Y[4]
    dY[2] = Y[5]
    r = np.sqrt(Y[0]**2+Y[1]**2+Y[2]**2)
    dY[3] = -mu*Y[0]/r**3
    dY[4] = -mu*Y[1]/r**3
    dY[5] = -mu*Y[2]/r**3
    return dY

def main():
    def derivFcn(t, y):
        return two_body_orbit(t, y, MU)

    r_0 = np.array([-5282.628, -4217.988, 296.511])
    v_0 = np.array([-4.541972, 5.885228, 2.043106])
    E = sf.cart2elm(r_0, v_0, MU)
    a = E[0]
    P = 2*m.pi*np.sqrt(a**3/MU)
    T0 = 0.0;
    times = np.arange(0, P*3, 20)
    tF = P*3
    dT = 20
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
    for i in range(np.size(r)):
        a[i] = np.divide(-MU, r[i]**2)
    t = output[:, 0]
    r_vec = np.empty([844, 3])
    v_vec = np.empty([844, 3])
    elms = np.empty([np.size(times)+1, 6])
    for i in range(np.size(r)):
        r_vec[i, 0] = output[i, 1]
        r_vec[i, 1] = output[i, 2]
        r_vec[i, 2] = output[i, 3]
        v_vec[i, 0] = output[i, 4]
        v_vec[i, 1] = output[i, 5]
        v_vec[i, 2] = output[i, 6]
        elms[i] = sf.cart2elm(r_vec[i, :], v_vec[i, :], MU)
    eng = .5*np.power(v, 2) - np.divide(MU, r) - (.5*lg.norm(v_0)**2 - MU/lg.norm(r_0))

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
    plt.xlabel("Hours", size=16)

    plt.figure()
    plt.title("Change in Energy vs. Time")
    plt.plot(t/3600, eng)
    plt.ylabel('Energy', size=16)
    plt.xlabel('Hours', size=16)

    plt.figure()
    plt.suptitle("Orbital Elements vs. Time")
    plt.subplot(3, 2, 1)
    plt.plot(t / 3600, elms[:, 0])
    a = plt.ylabel('a [km]')
    plt.subplot(3, 2, 2)
    plt.plot(t / 3600, elms[:, 1])
    b = plt.ylabel('e')
    plt.subplot(3, 2, 3)
    plt.plot(t / 3600, elms[:, 2])
    c = plt.ylabel('i [deg]')
    plt.subplot(3, 2, 4)
    plt.plot(t / 3600, elms[:, 3])
    d = plt.ylabel('RAAN [deg]')
    plt.subplot(3, 2, 5)
    plt.plot(t / 3600, elms[:, 4])
    e = plt.ylabel('Arg of Periapsis [deg]')
    plt.subplot(3, 2, 6)
    plt.plot(t / 3600, elms[:, 5])
    plt.xlabel('Hours')
    f = plt.ylabel('Anomaly [deg]')
    plt.show()
    return


if __name__ == "__main__":
    main()
