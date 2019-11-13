import numpy as np
import numpy.linalg as lg
import math as m
from scipy.integrate import ode


def orbit_prop(time_series, n, e, t_p):
    # allocate memory for anomalies
    E = np.empty(np.size(time_series))
    M = np.empty(np.size(time_series))
    theta = np.empty(np.size(time_series))
    e_norm = lg.norm(e)
    # Helper Functions
    def f(x, m_e):
        return x - e_norm * np.sin(x) - m_e
    def df(x):
        return 1 - e_norm * np.cos(x)
    # propagate through the time series
    for i in range(np.size(time_series)):
        M[i] = n*(time_series[i]-t_p)  # mean anomaly for this time step
        if M[i] < m.pi:  # inital guess based on mean anomaly
            guess = M[i]+e_norm/2
        else:
            guess = M[i]-e_norm/2
        it = 0
        error = 100.0
        while error > 10**-10 and it <= 50:  # newton raphson to find eccentric anomaly
            E[i] = guess-f(guess, M[i])/df(guess)
            error = np.abs((E[i]-guess)/E[i])
            guess = E[i]
            it = it+1
        theta[i] = 2*m.atan2(np.sqrt(1+e_norm)*np.tan(E[i]/2), np.sqrt(1-e_norm))  # find anomaly from eccentric anomaly
    return theta, E, M


def hyper_orbit_prop(time_series, n, e, t_p):
    # allocate memory for anomalies
    F = np.empty(np.size(time_series))
    M = np.empty(np.size(time_series))
    theta = np.empty(np.size(time_series))
    e_norm = lg.norm(e)
    # Helper Functions
    def f(x, m_h):
        return -x+e_norm*np.sinh(x)-m_h
    def df(x):
        return -1+e_norm*np.cosh(x)
    # propagate through the time series
    for i in range(np.size(time_series)):
        M[i] = n*(time_series[i]-t_p)  # mean anomaly for this time step
        if M[i] < m.pi:  # inital guess based on mean anomaly
            guess = M[i]+e_norm/2
        else:
            guess = M[i]-e_norm/2
        it = 0
        error = 100.0
        while error > 10**-10 and it <= 50:  # newton raphson to find eccentric anomaly
            F[i] = guess-f(guess, M[i])/df(guess)
            error = np.abs((F[i]-guess)/F[i])
            guess = F[i]
            it = it+1
        theta[i] = 2*m.atan2(np.sqrt(e_norm+1)*np.tanh(F[i]/2), np.sqrt(e_norm-1))  # find anomaly from eccentric anomaly
    return theta, F, M


def cart2elm(r, v, mu, deg=True):
    h = np.cross(r, v)
    r_norm = lg.norm(r)
    v_norm = lg.norm(v)
    e = np.cross(v, h) / mu - np.divide(r, r_norm)  # eccentricity
    e_norm = lg.norm(e)
    energy = (v_norm**2)/2 - mu/r_norm
    h_norm = lg.norm(h)
    k = (h_norm ** 2) / (r_norm * mu) - 1
    if energy < 0:
        a = -mu/(2*energy)
    elif -10e-12 < energy < 10e-12:
        a = m.inf
    else:
        a = mu/(2*energy)
    i = np.arccos(np.dot(h, [0, 0, 1])/h_norm)
    n = np.cross([0, 0, 1], h)
    n_norm = lg.norm(n)
    if e_norm < 10e-12 or e_norm > 10e-12:
        theta = np.arccos(k/e_norm)
        if np.dot(r,v)<0:
            theta = 2*m.pi-theta
        RAAN = np.arccos(np.dot(n, [1, 0, 0])/n_norm)
        omega = np.arccos(np.dot(n, e)/(e_norm*n_norm))
    if e_norm < 10e-12 and i < 10e-12:
        RAAN = 0
        omega = 0
        theta = np.arccos(r[1]/r_norm)
        if r[1] < 0:
            theta = 2*m.pi-theta
    elif e_norm < 10e-12:
        omega = 0
        RAAN = np.arccos(np.dot(n, [1, 0, 0]) / n_norm)
        theta = np.arccos(np.dot((n/n_norm),r)/r_norm)
        if r[2]< 0:
            theta = 2*m.pi-theta
    elif i < 10e-12:
        RAAN = 0
        omega = np.arccos(np.dot(e, [1, 0, 0])/e_norm)
        if e[1]< 0:
            omega = 2*m.pi-omega
    if deg:
        theta = 180*theta/m.pi
        i = 180*i/m.pi
        RAAN = 180*RAAN/m.pi
        omega = 180*omega/m.pi
    E = [a, e_norm, i, RAAN, omega, theta]
    return E


def elm2cart(E, mu, deg=True):
    # E - [a, e, i, RAAN, omega, theta]
    a = E[0]
    e = E[1]
    if deg:
        i = m.pi * E[2] / 180
        RAAN = m.pi * E[3] / 180
        omega = m.pi * E[4] / 180
        theta = m.pi * E[5] / 180
    p = a*(1 - e**2)
    r_pqw = np.array([(p/(1+e*np.cos(theta)))*np.cos(theta), (p/(1+e*np.cos(theta)))*np.sin(theta), 0])
    v_pqw = np.array([np.sqrt(mu/p)*(-np.sin(theta)), np.sqrt(mu/p)*(e+np.cos(theta)), 0])
    # R_3(-RAAN)R_1(-i)R_3(-omega)
    c1 = np.cos(-omega)
    c2 = np.cos(-i)
    c3 = np.cos(-RAAN)
    s1 = np.sin(-omega)
    s2 = np.sin(-i)
    s3 = np.sin(-RAAN)
    q1 = np.array([c1*c3-c2*s1*s3, c3*s1+c1*c2*s3, s3*s2])
    q2 = np.array([-c1*s3-c3*c2*s1, c1*c2*c3-s1*s3, c1*s2])
    q3 = np.array([s1*s2, -c1*s2, c2])
    Q = np.array([q1, q2, q3])
    r = np.matmul(Q, r_pqw)
    v = np.matmul(Q, v_pqw)
    return r, v

def R1(phi):
    return np.array([np.array([1, 0, 0 ]), np.array([0, np.cos(phi), np.sin(phi)]), np.array([0, -np.sin(phi), np.cos(phi)])])

def R2(phi):
    return np.array([np.array([np.cos(phi), 0, -np.sin(phi)]), np.array([0, 1, 0]), np.array([np.sin(phi), 0, np.cos(phi)])])

def R3(phi):
    return np.array([np.array([np.cos(phi), np.sin(phi), 0]), np.array([-np.sin(phi), np.cos(phi), 0]), np.array([0, 0, 1])])

def deg2rad(input):
    output = np.empty(np.size(input))
    for i in range(np.size(input)):
        output[i] = input[i]*np.pi/180
    return output

def orbit_prop_rk(r_0, v_0, T0, tF, dT):
    def two_body_orbit(t, Y, mu):
        dY = np.empty([6, 1])
        dY[0] = Y[3]
        dY[1] = Y[4]
        dY[2] = Y[5]
        r = np.sqrt(Y[0] ** 2 + Y[1] ** 2 + Y[2] ** 2)
        dY[3] = -mu * Y[0] / r ** 3
        dY[4] = -mu * Y[1] / r ** 3
        dY[5] = -mu * Y[2] / r ** 3
        return dY

    MU = 398600.4415
    RE = 6378.1363

    def derivFcn(t, y):
        return two_body_orbit(t, y, MU)

    Y_0 = np.concatenate([r_0, v_0], axis=0)
    rv = ode(derivFcn)

    #  The integrator type 'dopri5' is the same as MATLAB's ode45()!
    #  rtol and atol are the relative and absolute tolerances, respectively
    rv.set_integrator('dopri5', rtol=1e-10, atol=1e-20)
    rv.set_initial_value(Y_0, T0)
    output = []
    output.append(np.insert(Y_0, 0, T0))

    # Run the integrator and populate output array with positions and velocities
    while rv.successful() and rv.t < tF:  # rv.successful() and
        rv.integrate(rv.t + dT)
        output.append(np.insert(rv.y, 0, rv.t))

    #  Convert the output a numpy array for later use
    output = np.array(output)
    t = output[:, 0]

    r_vec = np.empty([np.shape(output)[0], 3])
    v_vec = np.empty([np.shape(output)[0], 3])

    for i in range(np.shape(output)[0]):
        r_vec[i, 0] = output[i, 1]
        r_vec[i, 1] = output[i, 2]
        r_vec[i, 2] = output[i, 3]
        v_vec[i, 0] = output[i, 4]
        v_vec[i, 1] = output[i, 5]
        v_vec[i, 2] = output[i, 6]
    return r_vec, v_vec
