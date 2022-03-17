# %% first clear all the variable that you


for name in dir():
    if not name.startswith('_'):
        del globals()[name]
del name

import sys
import os
sys.path.append("C:/Users/jiali/PANDA520-flowtube/PANDA520_flowtube/")
#file_path = "C:/Users/jiali/PANDA520-flowtube/PANDA520_flowtube/"
#os.chdir(file_path)

import numpy as np
import pandas as pd
from cmd_calib5 import cmd_calib5
from basic_input import inputs_va
from Vapour_calc import H2O_conc


class UnitFloat(float):

    def __new__(self, value, unit=None):
        return float.__new__(self, value)

    def __init__(self, value, unit=None):
        self.unit = unit


# summary for all the SA cali
# MION inlet
# 09.10 the frist SA cali with 41 cm for 3/4 inch and 58.5 cm for 1 inch first tower Br API9
# 10.28 the second SA cali with 50 cm for 3/4 inch and 68 cm for 1 inch we have Y pieces first tower Br API9
# 11.18 the third SA cali with 50 cm for 3/4 inch and 66 cm for 1 inch we have Y pieces and second tower Br API9
# 01.04 the fourth SA cali with 10 cm for 3/4 inch and 78 cm for 1 inch we have Y pieces and first tower NO3 tower API9
# NO3 inlet
# 01.27 the fifth SA cali with 26cm for 3/4 inch NO3 inlet API9
# MION inlet
# 02.07 the sixth SA cali with 10 cm for 3/4 inch and 61 cm for 1 inch first tower Br karsa
# %% Prepare the inputs
# T_cel = float(input('temperaure C:'))
# date = input('date:')
T_cel = 25
T = T_cel + 273.15  # K
p = 96060  # pressure Pa
date = ['09.10']
for i in range(len(date)):
    R1, L1, R2, L2, flag_tube, file, s1, s2 = inputs_va(date[i])

    # load H2O Q
    file = os.getcwd() + '/input_files/' + file

    H2O_data = pd.read_csv(file)

    if flag_tube == '3':
        Q1 = H2O_data['Q1'][0]  # lpm
        Q2 = H2O_data['Q2'][0]
        Q2 = Q1 + Q2
    elif flag_tube == '2':
        Q1 = H2O_data['Q1'][0]  # lpm
        Q2 = Q1
    else:
        Q1 = H2O_data['Q1'][0]
        Q2 = 0

    # It product for cali box
    Itx = 5.2009e10  # at Qx flow rate
    Qx = 20  # lpm

    N2Flow = 22  # slpm
    AirFlow = 50  # synthetic air slpm
    O2inAir = 0.209
    SO2Flow = 5  # slpm
    SO2BottlePpm = 5000  # ppm

    outflowLocation = 'before'  # outflow tube located before or after injecting air, water, and so2

    fullOrSimpleModel = 'full'  # simple: Gormley&Kennedy approximation, full: flow model (much slower)

    # calculate the conc for const species here, this file is SO2, O2, H2O
    WaterFlow1 = H2O_data['H2O_1']

    if outflowLocation in 'after':
        totFlow = N2Flow + AirFlow / 1000 + WaterFlow1 / 1000 + SO2Flow / 1000
    else:
        print(WaterFlow1.shape)
        totFlow = Q1 * np.ones(WaterFlow1.shape)

    O2conc1 = O2inAir * AirFlow / 1000 / totFlow * p / 1.3806488e-23 / T / 1e6
    H2Oconc1 = WaterFlow1 / 1000 / totFlow * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6
    SO2conc1 = SO2Flow / 1000 / totFlow * SO2BottlePpm * 1e-6 * p / 1.3806488e-23 / T / 1e6
    # H2Oconc1 = H2O_data['H2Oconc_1']
    WaterFlow2 = H2O_data['H2O_1']  # second H2O flow
    totFlow2 = Q2 * np.ones(WaterFlow1.shape)

    if flag_tube in ['1','2']:
        H2Oconc2 = H2Oconc1
        O2conc2 = O2conc1
        SO2conc2 = SO2conc2
    else:
        H2Oconc2 = (WaterFlow2) / 1000 / totFlow2 * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6
        # H2Oconc2 = H2O_data['H2Oconc_2']
        O2conc2 = O2conc1 * Q1 / Q2
        SO2conc2 = SO2conc1 * Q1 / Q2

    H2Oconc = np.transpose([H2Oconc1, H2Oconc2])
    O2conc = np.transpose([O2conc1, O2conc2])
    SO2conc = np.transpose([SO2conc1, SO2conc2])

    # % store all the const species to const_comp_conc
    const_comp_conc = np.transpose([SO2conc, H2Oconc, O2conc])

    if flag_tube == '3':
        const_comp_free = ['H2O', 'O2']
        const_comp_conc_free = [H2Oconc1[0], O2conc1[0]]
    else:
        const_comp_free = [0]
        const_comp_conc_free = [0]
    # calculate Initial concentration for some species here is OH and HO2
    csH2O = 7.22e-20  # cm2

    qyH2O = 1

    It = Itx * Qx / Q1

    OHconc = It * csH2O * qyH2O * H2Oconc1
    OHconc = OHconc.reset_index()['H2O_1']

    # store initial concentration
    Init_comp = ['OH', 'HO2']  # species have inital concentration
    Init_comp_conc = np.transpose([OHconc, OHconc])

    # diffusion constants, [cm**2/s]

    # DOH-air = 165 ± 20 Torr cm2 s-1, DHO2-He = 430 ± 30 Torr cm2 s-1，DO3-He = 410 ± 25 Torr cm2 s-1 at 296 K.
    # Source OH, HO2, and Ozone Gaseous Diffusion Coefficients
    dOH = 165 / (0.00750062 * p)  # convert to pa
    T0 = 296
    dOH = 101325 / p * dOH * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

    # Source A Measurement of the Diffusion Coefficient of Hydrogen Peroxide Vapor into Air
    dHO2 = 111 / (0.00750062 * p)  # convert to pa
    T0 = 296
    dHO2 = 101325 / p * dHO2 * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

    # https://www.engineeringtoolbox.com/air-diffusion-coefficient-gas-mixture-temperature-d_2010.html
    dH2O = 0.242  # cm2/s 20C
    T0 = 273.15 + 20
    dH2O = 101325 / p * dH2O * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

    # # can not find the source use the old one
    # dHSO3 = 0.126  # cm2/s
    # T0 = 300
    # dHSO3 = 101325 / p * dHSO3 * ((T ** (3 / 2)) / (T0 ** (3 / 2)))
    # % order is: HSO3, SO3, HO2, H2SO4, OH
    # D = [0.126 0.126 0.141 diff_sa_rh(298,rh)*1e4 0.215]

    # Diff_setname = ['HSO3','SO3','HO2','H2SO4','OH']
    # Diff_set = [0.126, 0.126, 0.141, 0.0895, 0.215]

    Diff_setname = ['OH', 'H2O', 'HO2', 'SO3']
    Diff_set = [dOH, dH2O, dHO2, 0.126]

    # grid for length and radius directions
    Zgrid = np.array(40).astype(int)  # number of grid points in tube length direction
    Rgrid = np.array(80).astype(int)  # number of grid points in tube radius direction

    # chemistry part
    sch_name = os.getcwd() + '/input_mechanism/SO2_SA.txt'
    chm_sch_mrk = ['{', 'RO2', '+', '', '', ';', '+', ';', '$', '{', ':', ';', '}']
    key_spe_for_plot = 'H2SO4'
    plot_spec = ['OH', 'HSO3', 'HO2', 'SO3', 'H2SO4']  # plot species
    const_comp = ['SO2', 'H2O', 'O2']  # species have constant concentration
    drh_str = str('0.*TEMP')
    erh_str = str('0.*TEMP')

    params = {'T': UnitFloat(T, "K"),  # temperature
              'p': UnitFloat(p, "Pa"),  # pressure pa
              'R1': UnitFloat(R1, "cm"),  # diameters for first tube
              'R2': UnitFloat(R2, "cm"),  # diameters for second tube
              'L1': UnitFloat(L1, "cm"),  # length for first tube
              'L2': UnitFloat(L2, "cm"),  # length for first tube
              'Q1': UnitFloat(Q1, "lpm"),  # flow for first tube
              'Q2': UnitFloat(Q2, "lpm"),  # flow for second tube
              'dt': 0.0001,  # flow for second tube
              'Diff_setname': Diff_setname,  # diffusion
              'Diff_set': Diff_set,
              # 'It': It, # it product for calculation
              'fullOrSimpleModel': fullOrSimpleModel,  # Gormley&Kennedy approximation, full: flow model (much slower)
              'sch_name': sch_name,  # file for the MCM file
              'chm_sch_mrk': chm_sch_mrk,  # markers to isolate sections of chemical scheme based on MCM KPP format
              'drh_str': str('0.*TEMP'),
              'erh_str': str('0.*TEMP'),
              'const_comp': const_comp,  # constant species
              'Init_comp': Init_comp,  # species have initial concentration
              'Zgrid': Zgrid,  # number of grid points in tube length direction
              'Rgrid': Rgrid,  # number of grid points in tube radius direction
              # 'formula': formula, # the formula for the plots
              'key_spe_for_plot': key_spe_for_plot,  # key species for plot
              'plot_spec': plot_spec,  # plot species
              'flag_tube': flag_tube,
              'const_comp_free': const_comp_free,
              'const_comp_conc_free': const_comp_conc_free
              }
    # %
    for i in range(8):
        print(list(params.keys())[i], list(params.values())[i], list(params.values())[i].unit)

    # i = 2
    # const_comp_conc= const_comp_conc[:,i,:]
    # Init_comp_conc=Init_comp_conc[i]
    # %% computation begins
    meanconc = []

    c = []

    for i in range(len(H2Oconc1)):  # range(H2SO4.size):
        if WaterFlow1[i] > 0:
            meanConc1, c1 = cmd_calib5(const_comp_conc[:, i, :], params, Init_comp_conc[i])
            meanconc.append(meanConc1)
            c.append(c1)

    meanconc_s = pd.DataFrame(np.transpose(meanconc))
    meanconc_s.index = plot_spec
    # %% save the modelled SA, HO2
    meanconc_s.to_csv('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_mean__' + s1)

    #with open('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_c' + s2, 'w') as f:
    #    for item in c:
    #        f.write("%s\n" % item)
    import csv
    with open('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_c' + s2, 'w') as f:
    # using csv.writer method from CSV package
         write = csv.writer(f)

         #write.writerow(fields)
         write.writerows(c)
# %% compare the data from matlab
# import scipy.io
# import matplotlib.pyplot as plt
# % order is: HSO3, SO3, HO2, H2SO4, OH
# data = scipy.io.loadmat(r'SA_calibrator_Miska/c_matlab.mat')
# c_mat = data['c']
# h2so4_m = c_mat[0][3]
# ho2_m = c_mat[0][2]

# h2so4_m = h2so4_m[:,-1]

# h2so4_p = c1[:,:,8][:,-1]


# plt.plot(h2so4_p,label = 'python')
# plt.plot(h2so4_m,label = 'matlab')
# plt.ylabel('conc')
# plt.xlabel('R')
# plt.legend()
# Ra = (h2so4_m-h2so4_p)
# plt.plot(Ra)
