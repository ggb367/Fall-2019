import numpy as np
import space_functions as sf
import numpy.linalg as lg

MU = 398600.4415
RE = 6378.1363
JDos = 0.00108248
f = 3.353e-3

pgd = 20.7984
lam = -156.3319
H = 3
thetaG = -23.6681 - lam
rho = 500
beta = 81.4
el = 25

[pgd, lam, thetaG, beta, el] = sf.deg2rad([pgd, lam, thetaG, beta, el])
r_ijk = [-6236.756, -844.821, 1789.052]
Q_ecf_enz = np.matmul(np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]), np.matmul(sf.R2(-pgd), sf.R3(lam)))

R_ecf = np.array([(RE/np.sqrt(1-(2*f-f**2)*np.sin(pgd)**2)+H)*np.cos(pgd)*np.cos(lam), (RE/np.sqrt(1-(2*f-f**2)*np.sin(pgd)**2)+H)*np.cos(pgd)*np.sin(lam), ((RE*(1-f)**2)/np.sqrt(1-(2*f-f**2)*np.sin(pgd)**2)+H)*np.sin(pgd)])
Q_ecf_ijk = sf.R3(thetaG)
rho_enz = np.matmul(Q_ecf_enz, np.matmul(Q_ecf_ijk, r_ijk)-R_ecf)
print(rho_enz)
print(lg.norm(rho_enz))
