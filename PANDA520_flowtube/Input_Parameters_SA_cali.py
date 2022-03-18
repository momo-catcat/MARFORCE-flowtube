# %% import package
# first clear all the variable that you
for name in dir():
    if not name.startswith('_'):
        del globals()[name]
del name

import sys

sys.path.append("C:/Users/jiali/PANDA520-flowtube/PANDA520_flowtube/")

import os
import numpy as np
import pandas as pd
from cmd_calib5 import cmd_calib5
from exp_setup import inputs_setup
import const_comp_conc_cal
import csv
from diffusion_const_added import add_diff_const as add_diff_const


# file_path = "C:/Users/jiali/PANDA520-flowtube/PANDA520_flowtube/"
# os.chdir(file_path)


# add unit after values
class UnitFloat(float):

    def __new__(self, value, unit=None):
        return float.__new__(self, value)

    def __init__(self, value, unit=None):
        self.unit = unit


'''
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
'''
# % Prepare the inputs
# T_cel = float(input('temperature C:'))
# date = input('date:')
T_cel = 25
T = T_cel + 273.15  # K
p = 96060  # pressure Pa
# select experiment name for simulating
date = ['09.10']  # can be multiple experiments
# %%
for i in range(len(date)):
    # %% load the parameters for experiment setup, change the basic_input
    R1, L1, R2, L2, flag_tube, file, s1, s2 = inputs_setup(date[i])

    # load H2O Q
    file = os.getcwd() + '/input_files/' + file
    H2O_data = pd.read_csv(file)
    # set the flow for tube 1 and tube 2 in each experiment
    if flag_tube in ['3', '4']:
        Q1 = H2O_data['Q1']  # lpm
        Q2 = H2O_data['Q2']
        Q2 = Q1 + Q2
    elif flag_tube == '2':
        Q1 = H2O_data['Q1']  # lpm
        Q2 = Q1
    else:
        Q1 = H2O_data['Q1']
        Q2 = Q1

    # set the constant precursors, change as you want
    outflowLocation = 'before'  # outflow tube located before or after injecting air, water, and so2

    fullOrSimpleModel = 'full'  # simple: Gormley&Kennedy approximation, full: flow model (much slower)
    # in this study, we have the outflow before injecting air, water and SO2

    const_comp_pre = ['SO2', 'O2']  # species always constant concentration
    const_comp_var = ['H2O']  # although constant for each exp but vary between different exps
    const_comp = const_comp_pre + const_comp_var
    N2Flow = 22  # lpm
    ppm = 1e-6

    # set the flows for the precursors
    pre_flow = [5, 50]  # slpm  AirFlow = 50  # synthetic air slpm  SO2Flow = 5  # slpm
    pre_standard_conc = [0.209, 5000 * ppm]  # O2inAir = 0.209 # SO2  SO2BottlePpm = 5000  # ppm

    # % store all the const species to const_comp_conc
    const_comp_conc, const_comp_free, const_comp_conc_free, OHconc = const_comp_conc_cal.const_comp_conc_cal(H2O_data,
                                                                                                             outflowLocation,
                                                                                                             const_comp_var,
                                                                                                             const_comp_pre, \
                                                                                                             pre_flow,
                                                                                                             N2Flow,
                                                                                                             pre_standard_conc, \
                                                                                                             Q1, Q2,
                                                                                                             T_cel, T,
                                                                                                             p,
                                                                                                             flag_tube)

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
              # 'Q1': UnitFloat(Q1, "lpm"),  # flow for first tube
              # 'Q2': UnitFloat(Q2, "lpm"),  # flow for second tube
              'dt': 0.0001,  # flow for second tube
              'Diff_setname': Diff_setname,  # diffusion for the species that you want to have
              'Diff_set': Diff_set,
              # 'It': It, # it product for calculation
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
    for i in range(6):
        print(list(params.keys())[i], list(params.values())[i], list(params.values())[i].unit)

    # %% computation begins
    meanconc = []

    c = []

    for i in range(len(H2O_data)):
        if H2O_data['H2O_1'][i] > 0:
            meanConc1, c1 = cmd_calib5(const_comp_conc[:, :, i], params, Init_comp_conc[i], Q1[i], Q2[i])
            meanconc.append(meanConc1)
            c.append(c1)

    meanconc_s = pd.DataFrame(np.transpose(meanconc))
    meanconc_s.index = plot_spec
    # %% save the modelled SA, HO2
    meanconc_s.to_csv('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_mean__' + s1)

    with open('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_c' + s2, 'w') as f:
        # using csv.writer method from CSV package
        write = csv.writer(f)

        write.writerows(c)
