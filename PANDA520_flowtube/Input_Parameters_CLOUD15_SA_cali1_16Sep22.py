# this file is used for all the inputs

# clear variables
for name in dir():
    if not name.startswith('_'):
        del globals()[name]
del name

# %% import packages and functions
import sys
# sys.path.append('/Users/hexuchen/Documents/science/coding/git/Projects/PANDA520-flowtube/PANDA520_flowtube/')
import os
import numpy as np
import pandas as pd
import csv
from cmd_calib5 import cmd_calib5
from exp_setup import inputs_setup_CLOUD15
from Calcu_by_flow import const_comp_conc_cal, const_comp_conc_cal_H2O, const_comp_conc_cal_OH
from diffusion_const_added import add_diff_const as add_diff_const
from Vapour_calc import H2O_conc


# file_path = '/Users/hexuchen/Documents/science/coding/git/Projects/PANDA520-flowtube/PANDA520_flowtube/'
# os.chdir(file_path)

# add unit after values
class UnitFloat(float):

    def __new__(self, value, unit=None):
        return float.__new__(self, value)

    def __init__(self, value, unit=None):
        self.unit = unit


'''
# summary for all the SA cali
CLOUD15 SA calibration using UH calibrator 1.
'''

# set temperature and pressure
T_cel = 30.5
T = T_cel + 273.15  # K
p = 101000  # pressure Pa

# select experiment name for simulating
date = ['CLOUD15_cali1_calibrator1_16Sep22']  # can be multiple experiments

#%%
for i in range(len(date)):

    # %% load the parameters for experiment setup, change the basic_input
    R1, L1, R2, L2, flag_tube, file, s1, s2, H2O_1, H2O_2, H2Oconc_1, H2Oconc_2, Q1, Q2 = inputs_setup_CLOUD15(date[i])


    # if H2Oconc_1[0] == 1: #if there is no H2O concentration calculate based on the flows. Temp setting
    #     H2Oconc_1 = H2O_conc(T_cel, H2O_1 / 1000 / (H2O_1 / 1000 + Q1))
    #     H2Oconc_1 = H2Oconc_1.H2O_conc.values[:]
    #     H2Oconc_2 = H2O_conc(T_cel, H2O_2 / 1000 / (H2O_2 / 1000 + Q2))
    #     H2Oconc_2 = H2Oconc_2.H2O_conc.values[:]


    # set the constant precursors, change as you want
    outflowLocation = 'before'  # outflow tube located before or after injecting air, water, and so2

    fullOrSimpleModel = 'full'  # simple: Gormley&Kennedy approximation, full: flow model (much slower)
    # in this study, we have the outflow before injecting air, water and SO2
    O2flow = 100.0  # slpm
    SO2flow = 10.7  # slpm
    sampflow = 20.0  # lpm
    ppm = 1e-6
    O2ratio = 0.209  # O2inAir = 0.209
    SO2ratio = 1000 * ppm
    N2Flow = 24.0  # lpm
    const_comp_pre = ['SO2', 'O2']  # species have constant concentration and are calculated from flows
    const_comp_pre_know = ['H2O']  # species have known constant concentration but already known
    const_comp = const_comp_pre + const_comp_pre_know  # species have constant concentration
    # get all the concentrations
    O2conc, SO2conc = const_comp_conc_cal(O2flow, SO2flow,outflowLocation, sampflow, H2O_1,  N2Flow, O2ratio, SO2ratio,
                                 Q1, Q2, T_cel, T, p, flag_tube)


    H2Oconc = const_comp_conc_cal_H2O(O2flow, SO2flow, outflowLocation, sampflow, H2O_1, H2O_2, N2Flow, O2ratio,
                                      Q1, Q2, T_cel, T, p, flag_tube)

    print('H2O concentrations')
    print(H2Oconc)

    # % store all the const species to const_comp_conc follow the order of const_comp
    const_comp_conc = np.transpose([SO2conc, O2conc, H2Oconc])

    OHconc, const_comp_free, const_comp_conc_free = const_comp_conc_cal_OH(H2Oconc, O2conc, Q1, flag_tube)

    # store initial concentration
    Init_comp = ['OH', 'HO2']  # species have inital concentration
    Init_comp_conc = np.transpose([OHconc, OHconc])

    # add some diffusion constants add more in diffusion_const_added.py file if you want
    Diff_setname = ['OH', 'H2O', 'HO2', 'SO3']
    dOH, dH2O, dHO2 = add_diff_const(p, T)
    Diff_set = [dOH, dH2O, dHO2, 0.126]

    # grid for length and radius directions
    Zgrid = np.array(40).astype(int)  # number of grid points in tube length direction
    Rgrid = np.array(80).astype(int)  # number of grid points in tube radius direction

    # chemistry part
    sch_name = os.getcwd() + '/input_mechanism/SO2_SA.txt'
    chm_sch_mrk = ['{', 'RO2', '+', '', '', ';', '+', ';', '$', '{', ':', ';', '}']
    key_spe_for_plot = 'H2SO4'
    plot_spec = ['OH', 'HSO3', 'HO2', 'SO3', 'H2SO4']  # plot species

    params = {'T': UnitFloat(T, "K"),  # temperature
              'p': UnitFloat(p, "Pa"),  # pressure pa
              'R1': UnitFloat(R1, "cm"),  # diameters for first tube
              'R2': UnitFloat(R2, "cm"),  # diameters for second tube
              'L1': UnitFloat(L1, "cm"),  # length for first tube
              'L2': UnitFloat(L2, "cm"),  # length for first tube
              'dt': 0.0001,  # dt * timesteps * numLoop is time elapsed in the final solution
              'Diff_setname': Diff_setname,  # diffusion for the species that you want to have
              'Diff_set': Diff_set,
              'fullOrSimpleModel': fullOrSimpleModel,  # Gormley&Kennedy approximation, full: flow model (much slower)
              'sch_name': sch_name,  # file for the MCM file
              'chm_sch_mrk': chm_sch_mrk,  # markers to isolate sections of chemical scheme based on MCM KPP format
              'const_comp': const_comp,  # constant species
              'Init_comp': Init_comp,  # species have initial concentration
              'Zgrid': Zgrid,  # number of grid points in tube length direction
              'Rgrid': Rgrid,  # number of grid points in tube radius direction
              # 'formula': formula, # the formula for the plots
              'key_spe_for_plot': key_spe_for_plot,  # key species for plot and criterion to stop the loop
              'plot_spec': plot_spec,  # plot species
              'flag_tube': flag_tube,
              'const_comp_free': const_comp_free,
              'const_comp_conc_free': const_comp_conc_free
              }
    # %
    # for j in range(6):
    #     print(list(params.keys())[i], list(params.values())[i], list(params.values())[i].unit)

    # %% computation begins
    meanconc = []

    c = []

    for j in range(len(H2O_1)):
        if OHconc[j] > 0:
            meanConc1, c1 = cmd_calib5(const_comp_conc[:, j, :], params, Init_comp_conc[j], Q1[j], Q2[j])
            meanconc.append(meanConc1)
            print(meanconc)
            c.append(c1)

    meanconc_s = pd.DataFrame(meanconc,columns=plot_spec)
    # meanconc_s.index = plot_spec


    # % save the modelled SA, HO2
    meanconc_s.to_csv(os.getcwd() + '/Export_files/SA_cali_' + str(date[0]) + '.csv')

    with open(os.getcwd() + '/Export_files/SA_cali_' + str(date[0]) + '.txt', 'w') as f:
        # using csv.writer method from CSV package
        write = csv.writer(f)

        write.writerows(c)
