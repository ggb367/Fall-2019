import numpy as np
import space_functions as sf

MU = 398600.4415
RE = 6378.1363
JDos = 0.00108248
f = 3.353e-3

pgd = 20.7984
lam = -156.3319
H = 3000
thetaG = -23.6681 - lam
rho = 500
beta = 81.4
el = 25

[pgd, lam, thetaG, beta, el] = sf.deg2rad([pgd, lam, thetaG, beta, el])

Q_ecf_enz = np.matmul(np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]), np.matmul(sf.R2(-pgd), sf.R3(lam)))
rho_enz = np.array([rho*np.cos(el)*np.sin(beta), rho*np.cos(el)*np.cos(beta), rho*np.sin(el)])
R_ecf = np.array([(RE/np.sqrt(1-(2*f-f**2)*np.sin(pgd)**2)+H)*np.cos(pgd)*np.cos(lam), (RE/np.sqrt(1-(2*f-f**2)*np.sin(pgd)**2)+H)*np.cos(pgd)*np.sin(lam), ((RE*(1-f)**2)/np.sqrt(1-(2*f-f**2)*np.sin(pgd)**2)+H)*np.sin(pgd)])
Q_ecf_ijk = sf.R3(thetaG)
r_ijk = np.matmul(np.transpose(Q_ecf_ijk), (np.matmul(np.transpose(Q_ecf_enz), rho_enz)+R_ecf))
print(r_ijk)

