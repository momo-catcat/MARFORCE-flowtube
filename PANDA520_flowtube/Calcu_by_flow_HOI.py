# this file is used for HOI calibration with Y piece for dilution

import numpy as np
import pandas as pd
import os
import sys
from pandas import Index
from Vapour_calc import H2O_conc

#%%
def calculate_concs(paras):
    # %%
    dirpath = os.path.dirname(__file__)  # get the current path of this file
    #dirpath = '/Users/jiali/Documents/work/AMT-MION2-Paper/PANDA520_flowtube' # for test

    folder_flowtube = dirpath + '/'
    sys.path.append(folder_flowtube)
    ## Input file folder
    input_file_folder = dirpath + '/input_files/'
    input_mechanism_folder = dirpath + '/input_mechanism/'
    ## Add the folder where the final results will be exported
    export_file_folder = dirpath + '/export_files/'
    file_name = paras['file_name']
    data = pd.read_csv(input_file_folder + file_name)
    inputs: Index = data.columns.tolist()
    outflowLocation = paras['outflowLocation']
    sampflow = paras['sampleflow']
    O2ratio = paras['O2ratio']
    p = paras['p']
    T = paras['T']
    Itx = paras['Itx']
    Qx = paras['Qx']

    # %%
    flag_tube12 = ['H2Oflow',
                   'N2flow']  # flag_tube '1' and '2' should include these two columns otherwise, flag_tube3 = ['H2Oflow1','H2Oflow2','N2flow1','N2flow2']
    # check all the flows
    if all(x in inputs for x in flag_tube12):
        N2flow = data['N2flow']
        O2flow = data['O2flow']
        I2conc = data['I2conc']
        I2flow = data['I2flow']
        H2Oflow = data['H2Oflow']
        Q = data['Q']
        sumflow = O2flow + I2flow + H2Oflow + N2flow
        It = Itx * Qx / Q * 1e3
        if 'L2' in paras.keys():
            flag_tube = '2'
        else:
            flag_tube = '1'
            paras['L2'] = np.array(0, dtype=np.float64)
            paras['R2'] = paras['R1']
    else:
        O2flow1 = data['O2flow1']
        O2flow2 = data['O2flow2']
        H2Oflow1 = data['H2Oflow1']
        H2Oflow2 = data['H2Oflow2']
        N2flow1 = data['N2flow1']
        N2flow2 = data['N2flow2']
        I2conc1 = data['I2conc1']
        I2conc2 = data['I2conc2']
        Q1 = data['Q1']
        Q2 = data['Q2']
        sumflow1 = O2flow1 +  H2Oflow1 + N2flow1
        It = Itx * Qx / Q1 *1e3
        flag_tube = '3'
#%%
    # check if we have the H2O concentration
    if all(x in inputs for x in ['H2Oconc']):
        H2O_concentration = data['H2Oconc']
    elif all(x in inputs for x in ['H2Oconc1', 'H2Oconc2']):
        H2Oconc_1 = data['H2Oconc1']
        H2Oconc_2 = data['H2Oconc2']

    kB = 1.3806488e-23  # boltzmann constant

    if outflowLocation in 'after':
        totFlow1 = sumflow
    elif flag_tube in ['3']:
        totFlow1 = Q1
        totFlow2 = sampflow *1e3
        H2Oconc1 = H2Oflow1  / totFlow1 * H2O_conc(T, 1).SatP[0] / kB / T / 1e6
        O2conc1 = O2flow1 * O2ratio  / totFlow1 * p / kB / T / 1e6
    else:
        totFlow1 = sampflow * 1e3
        H2Oconc1 = H2Oflow/ totFlow1 * H2O_conc(T, 1).SatP[0] / kB / T / 1e6
        O2conc1 = O2flow * O2ratio / totFlow1 * p / kB / T / 1e6



    if flag_tube in ['1', '2']:
        O2conc2 = O2conc1

        H2Oconc2 = H2Oconc1

        I2conc1 = I2conc
        I2conc2 = I2conc
    else:
        O2conc2 = (O2flow1 + O2flow2) * O2ratio  / totFlow2 * p / kB / T / 1e6
        H2Oconc2 = (H2Oflow1 + H2Oflow2)  / totFlow2 * H2O_conc(T, 1).SatP[0] / kB / T / 1e6

    csH2O = 7.22e-20  # cm2
    qyH2O = 1

    O2conc = np.transpose([O2conc1, O2conc2])
    I2conc = np.transpose([I2conc1, I2conc2])
    #%%
    if 'H2Oconc_1' in locals():
        H2Oconc = np.transpose([H2Oconc_1, H2Oconc_2])
        OHconc = It * csH2O * qyH2O * H2Oconc_1
    elif 'H2O_concentration ' in locals():
        H2Oconc = np.transpose([H2O_concentration , H2O_concentration])
        OHconc = It * csH2O * qyH2O * H2O_concentration
    else:
        H2Oconc = np.transpose([H2Oconc1, H2Oconc2])
        OHconc = It * csH2O * qyH2O * H2Oconc1

    if flag_tube == '3':
        const_comp_free = ['H2O', 'O2']
        const_comp_conc_free = [H2Oconc[:, 0], O2conc[:, 0]]
    else:
        const_comp_free = []
        const_comp_conc_free = [0]


    paras['Q1'] = Q1
    paras['Q2'] = Q2

    paras['const_comp_free'] = const_comp_free
    paras['const_comp_conc_free'] = const_comp_conc_free
    paras['sch_name'] = input_mechanism_folder + paras['sch_name']
    paras['OHconc'] = OHconc
    paras['flag_tube'] = flag_tube

    #if 'N2flow' in locals():
    #    paras['N2flow'] = N2flow
    #else:
    #    paras['N2flow1'] = N2flow1
    #    paras['N2flow2'] = N2flow2
    print('Q1',Q1)
    print('Q2',Q2)
    # %%
    return O2conc, I2conc, H2Oconc, paras, export_file_folder
