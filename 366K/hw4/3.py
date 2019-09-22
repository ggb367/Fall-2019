import numpy as np
import numpy.linalg as lg
import math
import matplotlib.pyplot as plt

# constants
MU = 398600.4415  # DU^3/TU^2

def orbit_prop(time_series, n, e, t_p):
    MAX_ITER = 50
    E = np.empty(np.size(time_series))
    M = np.empty(np.size(time_series))
    θ = np.empty(np.size(time_series))
    for i in range(np.size(time_series)):
        m_e = n*(time_series[i]-t_p)
        e_norm = lg.norm(e)

        def f(x):
            return x-e_norm*np.sin(x)-m_e

        def df(x):
            return 1-e_norm*np.cos(x)
        if m_e < math.pi:
            guess = m_e+e_norm/2
        else:
            guess = m_e-e_norm/2
        tol = 0.0000000001
        it = 0
        error = 100.0
        while error > tol and it <= MAX_ITER:
            thisE = guess - f(guess) / df(guess)
            error = np.abs((thisE - guess)/thisE)
            guess = thisE
            it = it+1
        E[i] = thisE
        M[i] = m_e
        θ[i] = 2*math.atan2(np.sqrt(1+e_norm)*np.tan(thisE/2), np.sqrt(1-e_norm))
    return θ, E, M

# Initial conditions
r_0 = [-17310, 39, -784]  # km
v_0 = [-.251, -2.827, 3.754]  # km/s

# a
r_0_norm = lg.norm(r_0)
v_0_norm = lg.norm(v_0)
h = np.cross(r_0, v_0)
e = np.cross(v_0, h)/MU - np.divide(r_0, r_0_norm)
e_norm = lg.norm(e)
if e_norm == 0:
    print('a. It is a circular orbit')
    print("e is "+str(e_norm))
elif e_norm < 1:
    print("a. It is an elliptical orbit")
    print("e is "+str(e_norm))
elif e_norm == 1:
    print("a. It is a parabolic orbit")
    print("e is "+str(e_norm))
elif e_norm > 1:
    print("a. It is a hyperbolic orbit")
    print("e is "+str(e_norm))

# b
h_norm = lg.norm(h)
k = (h_norm*h_norm)/(r_0_norm*MU) - 1
theta_0 = np.arccos(k/e_norm)
if np.dot(r_0, v_0) < 0:
    theta_0 = 2*np.pi() - theta_0
print("b. θ at t_0 is "+str(theta_0))
E_0 = 2*math.atan2(np.sqrt(1-e_norm)*np.tan(theta_0/2), np.sqrt((1+e_norm)))
print("Eccentric Anomaly at t_0 is: "+str(E_0))
M_E_0 = E_0-e_norm*np.sin(E_0)
print("Mean Anomaly at t_0 is: "+str(M_E_0))
# c
a = (r_0_norm*(1+e_norm*np.cos(theta_0)))/(1-e_norm*e_norm)
n = np.sqrt(MU/math.pow(a, 3))
print(h_norm)
t_p = M_E_0/n
print("c. t_p is: "+str(t_p))

# d
t = np.linspace(0, 40000, 81)
[theta, E, M] = orbit_prop(t, n, e, t_p)
cos = np.cos(theta)
r = np.multiply((h_norm*h_norm/MU), np.divide(1, (1+np.multiply(e_norm, cos))))  # orbit equation
a = (r[8]*(1+e_norm*np.cos(theta[8])))/(1-e_norm*e_norm)  # use any r and theta paring to find a
v = np.sqrt(MU*(np.divide(2, r))-1/a)  # use r and a to find v, from vis-viva
fig, axs = plt.subplots(3, sharex=True, gridspec_kw={'hspace': .3})
axs[0].plot(t, theta)
axs[0].set_title("Anomaly vs Time")
axs[1].plot(t, r)
axs[1].set_title("Radius vs Time")
axs[2].plot(t, v)
axs[2].set_title("Velocity vs Time")
plt.xlabel('Time [sec]')
fig.suptitle('Orbit Propagation', fontsize=16, fontweight='bold')
plt.show()
