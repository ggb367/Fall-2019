import numpy as np
import numpy.linalg as lg
import space_functions as sf

MU = 398600.4415
RE = 6378.1363
JDos = 0.00108248
f = 3.353e-3

E_i = [8000, 0.1, 50, 45, 300, 0]
E_f = [8500, 0.15294, 50, 45, 300, 0]
r_i, v_i = sf.elm2cart(E_i, MU)
r_f, v_f = sf.elm2cart(E_f, MU)
delta_v = v_f-v_i
print("∆v:")
print(delta_v)
v_r = np.dot(delta_v, r_i/lg.norm(r_i))
v_p = np.sqrt(lg.norm(delta_v)**2-v_r**2)
print("|∆v_r|: "+str(v_r)+" km/s |∆v_⊥|: "+str(v_p)+" km/s")
