import numpy as np
import numpy.linalg as lg
import math
import matplotlib.pyplot as plt

# constants
MU = 398600.4415  # DU^3/TU^2

def hyper_orbit_prop(time_series, n, e, t_p):
    # allocate memory for anomalies
    F = np.empty(np.size(time_series))
    M = np.empty(np.size(time_series))
    θ = np.empty(np.size(time_series))
    e_norm = lg.norm(e)
    # Helper Functions
    def f(x, m_h):
        return -x+e_norm*np.sinh(x)-m_h
    def df(x):
        return -1+e_norm*np.cosh(x)
    # propagate through the time series
    for i in range(np.size(time_series)):
        M[i] = n*(time_series[i]-t_p)  # mean anomaly for this time step
        if M[i] < math.pi:  # inital guess based on mean anomaly
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
        θ[i] = 2*math.atan2(np.sqrt(e_norm+1)*np.tanh(F[i]/2), np.sqrt(e_norm-1))  # find anomaly from eccentric anomaly
    return θ, F, M

# Initial conditions
r_0 = [0, 39991, -52979]  # km
v_0 = [0, 0, 3.5511]  # km/s

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
theta_0 = np.arccos(k/e_norm)  # anomaly at initial condition
if np.dot(r_0, v_0) < 0:
    theta_0 = 2*math.pi - theta_0
print("b. Anomaly at t_0 is "+str(theta_0))
F_0 = 2*math.atanh(np.sqrt((e_norm-1)/(e_norm+1))*np.tan(theta_0/2))
print("Eccentric Anomaly at t_0 is: "+str(F_0))
M_h_0 = -F_0+e_norm*np.sinh(F_0)
print("Mean Anomaly at t_0 is: "+str(M_h_0))

# c
a = h_norm**2/MU*(e_norm**2-1)**-1
n = np.sqrt(MU/a**3)
t_p = -M_h_0/n
print("c. t_p is: "+str(t_p))

# d
t = np.arange(-20000, 20000, 500)
t_hour = np.divide(t, 3600)
theta, E, M = hyper_orbit_prop(t, n, e, t_p)
print("Anomaly at t=10000 is: "+str(theta[60]))
r = np.multiply((h_norm*h_norm/MU), np.divide(1, (1+np.multiply(e_norm, np.cos(theta)))))  # orbit equation
a = (r[8]*(1+e_norm*np.cos(theta[8])))/(1-e_norm*e_norm)  # use any r and theta paring to find a
v = np.sqrt(MU*(np.divide(2, r))-1/a)  # use r and a to find v, from vis-viva
KE = 0.5*np.square(v)
PE = -np.divide(MU, r)
TE = KE+PE
dTE = TE - (KE[0]+PE[0])
fig, axs = plt.subplots(3, sharex=True, gridspec_kw={'hspace': .3})
axs[0].plot(t_hour, theta)
axs[0].set_title("Anomaly vs Time")
axs[0].set_ylabel("Anomaly [rad]")
axs[1].plot(t_hour, r)
axs[1].set_title("Radius vs Time")
axs[1].set_ylabel("Radius [km]")
axs[2].plot(t_hour, v)
axs[2].set_title("Velocity vs Time")
axs[2].set_ylabel("Velocity [km/s]")
plt.xlabel('Time [Hours]')
fig.suptitle('Orbit Propagation', fontsize=16, fontweight='bold')
plt.figure(2)
plt.plot(t_hour, dTE)
plt.title("Change in Energy vs Time")
plt.xlabel('Time [Hours]')
plt.ylabel('Specific Energy [kJ/kg]')
plt.ylim([-0.5, 0.5])
plt.show()
