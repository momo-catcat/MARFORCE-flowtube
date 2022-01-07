# This function is used to calculate theoretical concentration profile of particles by Alonso et al. 2010 (e.g. 23)
# Ref: 10.1016/j.jaerosci.2010.02.003

import numpy as np
from pene_rate import cal_pene_rate as cal_pene_rate
import matplotlib.pyplot as plt

def conc_profile_Alonso(R, r, MM, rho, L, Q, P = 101000, T = 298.15):
    #R -- radius of the tube (m)
    #r -- radius of grads (m)
    #MM -- molecular mass of given molecule (g mol-1)
    #rho -- bulk density of given material (kg m-3)
    #L -- length of tube (m)
    #Q -- flow rate (m3 s-1)
    #x -- calculate the penetration rate at axial position x (m)

    pene_rate = cal_pene_rate(P, T, MM, rho, L, Q, R)

    D = pene_rate.D[0]


    D = 0.0921 * 1e-4 #adopted from matlab

    u = Q / (np.pi * R ** 2)

    norm_r = r / R

    beta = D * L / u / R ** 2

    print(beta)

    # for beta * x >= 0.2
    conc_prof = 1.477 * np.exp(-3.659 * beta) * (1 - 1 / 11 * (18 * norm_r ** 2 - 9 * norm_r ** 4 + 2 * norm_r
                                                                   ** 6))

    return(conc_prof)


R = 0.0125
conc_prof = conc_profile_Alonso(R, np.linspace(0, R, 100), 176, 4600, 2, 22 / 60 / 1000)

# print(conc_prof)

plt.plot(np.linspace(0, R, 100), conc_prof)
plt.show()