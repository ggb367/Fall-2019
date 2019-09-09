import numpy as np
import matplotlib.pyplot as plt

C = 172800/(60*60*24) #sec

t = np.linspace(0, 10, 10000)
f = np.multiply(2*np.pi, t) #2*np.pi()*t
f = np.cos(np.divide(f, C))
plt.plot(t, f)
plt.xlabel('Time(days)')
plt.ylabel('f(t)')
plt.title('f(t) = cos(2*pi*t/C)')
plt.show()