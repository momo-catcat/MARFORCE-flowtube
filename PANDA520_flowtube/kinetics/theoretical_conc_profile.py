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

    print(D)


    u = Q / (np.pi * R ** 2)

    norm_r = r / R

    beta = D * L / u / R ** 2

    print('Beta = ' + str(beta) + '. Note Beta needs to be larger than 0.2')

    # for beta * x >= 0.2
    conc_prof = 1.477 * np.exp(-3.659 * beta) * (1 - 1 / 11 * (18 * norm_r ** 2 - 9 * norm_r ** 4 + 2 * norm_r ** 6))

    return(conc_prof)


R = 0.005
conc_prof = conc_profile_Alonso(R, np.linspace(0, R, 100), 97, 1830, 2, 8.5  / 60 / 1000)
## import theoretical value
import pandas as pd
prof_conc_theory = pd.read_csv('../Export_files/Theoretical_model.csv')

plt.plot(np.linspace(0, R, 100), conc_prof)
plt.plot(prof_conc_theory.R/100, prof_conc_theory.SA  / 1e8)
plt.legend(['Theory','Model'])
plt.show()


# ##data export
# file_path = '/Users/momo/Documents/science/publication/IO3_measurement/data/Figures/Model_kinetic_validation.xlsx'
# Export_data = pd.concat([pd.DataFrame({'Theory_R': np.linspace(0,R,100)}),
#                            pd.DataFrame({'Theory_SA_remain': conc_prof}),
#                            pd.DataFrame({'Model_R': prof_conc_theory.R/100}),
#                            pd.DataFrame({'Model_SA_remain': prof_conc_theory.SA / 1e8})], axis = 1)
#
# Export_data.to_excel(file_path)