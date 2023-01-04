import sys
import pandas as pd
import numpy as np
import os



##--------/* Add path */--------
## Add the folder where the software package locates (The PANDA520_flowtube folder)
dirpath = os.path.dirname(__file__) #get the current path of this file
#dirpath = 'C:/Users/jiali/PANDA520-flowtube/PANDA520_flowtube'
folder_flowtube = dirpath + '/'
sys.path.append(folder_flowtube)
## Input file folder
input_file_folder = dirpath + '/input_files_old/'
## Add the folder where the final results will be exported
export_file_folder= dirpath + '/export_files/'

#%%
##--------/* Information of tubes and water content & flows for different calibration stages */--------
def inputs_setup(date):
    # flag_tube, '4','3','2','1' refers to the setup of the experiment
    # '4' two tubes with different inner diameters and have Y piece, run the second tube with two flows simultaneously
    # '3' same as '4', but run the second tube with one flow after converting the mean concentrations
    # '2' two tubes with different inner diameters
    # '1' one tube
    # R1, inner diameter for the first tube, unit cm
    # L1, length for the first tube, unit cm
    # R2, inner diameter for the second tube if there is any, = 0 for '1', unit cm
    # L2, length for the second tube, unit cm
    # Itx, It product at Qx flow rate
    # SO2flow, SO2 flow, sccm
    # O2flow, synthetic air flow, sccm
    # N2flow, main N2 flow, lpm
    # H2O_1, H2O flow for the first tube, sccm
    # H2O_2, H2O flow for the second tube, sccm
    # T, temperature, Celsius degree
    # Q1, flow for the first tube, unit lpm
    # Q2, extra flow for the second tube (Q1+Q2 is the flow for second tube), unit lpm
    # file, the file containing information of water, UV and Q, for example, 'H2O_1.csv'
    if date == 'CLOUD15_cali1_calibrator1_16Sep22':
            flag_tube = '2'
            R1 = 0.78
            L1 = 13
            R2 = 1.04
            L2 = 60
            Itx = 5.2009e10
            Qx = 20
            file = 'H2O_CLOUD15_16Sep22.csv'
            file = input_file_folder + file
            H2O_data = pd.read_csv(file, delimiter = ';')
            H2O_1 = H2O_data['H2O_1']
            H2O_2 = H2O_data['H2O_2']
            H2Oconc_1 = H2O_data['H2Oconc_1']
            H2Oconc_2 = H2O_data['H2Oconc_2']
            SO2flow = H2O_data['SO2flow']
            O2flow = H2O_data['O2flow']
            N2flow = H2O_data['N2flow']
            T_cel = H2O_data['T']
            Q1 = H2O_data['Q1']
            Q2 = H2O_data['Q1'] + H2O_data['Q2']
    ## We can process multiple days at the same time, just follow the same codes as above
    ## For example, we use elif as below:
    elif date == 'CLOUD15_cali2_calibrator2_18Sep22':
            flag_tube = '2'
            R1 = 0.5
            L1 = 17
            R2 = 1.04
            L2 = 60
            Itx = 6.186e11*0.1
            Qx = 7.6
            file = 'H2O_CLOUD15_18Sep22.csv'
            file = input_file_folder + file
            H2O_data = pd.read_csv(file, delimiter = ',')
            H2O_1 = H2O_data['H2O_1']
            H2O_2 = H2O_data['H2O_2']
            H2Oconc_1 = H2O_data['H2Oconc_1']
            H2Oconc_2 = H2O_data['H2Oconc_2']
            SO2flow = H2O_data['SO2flow']
            O2flow = H2O_data['O2flow']
            N2flow = H2O_data['N2flow']
            T_cel = H2O_data['T']
            Q1 = H2O_data['Q1']
            Q2 = H2O_data['Q1'] + H2O_data['Q2']
    elif date == 'CLOUD15_Br_UHL_04Nov22':
            flag_tube = '2'
            R1 = 0.78
            L1 = 13
            R2 = 1.04
            L2 = 60
            Itx = 5.2009e10 #6.186e11*0.1
            Qx = 20
            file = 'H2O_CLOUD15_04Nov22_Br.csv'
            file = input_file_folder + file
            H2O_data = pd.read_csv(file, delimiter = ',')
            H2O_1 = H2O_data['H2O_1']
            H2O_2 = H2O_data['H2O_2']
            H2Oconc_1 = H2O_data['H2Oconc_1']
            H2Oconc_2 = H2O_data['H2Oconc_2']
            SO2flow = H2O_data['SO2flow']
            O2flow = H2O_data['O2flow']
            N2flow = H2O_data['N2flow']
            T_cel = H2O_data['T']
            Q1 = H2O_data['Q1']
            Q2 = H2O_data['Q1'] + H2O_data['Q2']
    elif date == 'CLOUD15_NO3_UHL_03Nov22':
            flag_tube = '2'
            R1 = 0.78
            L1 = 10.5
            R2 = 0.5
            L2 = 28.5
            Itx = 5.2009e10 # 6.186e11 * 0.1
            Qx = 20
            file = 'H2O_CLOUD15_04Nov22_NO3_new_attenuator.csv'
            file = input_file_folder + file
            H2O_data = pd.read_csv(file, delimiter=',')
            H2O_1 = H2O_data['H2O_1']
            H2O_2 = H2O_data['H2O_2']
            H2Oconc_1 = H2O_data['H2Oconc_1']
            H2Oconc_2 = H2O_data['H2Oconc_2']
            SO2flow = H2O_data['SO2flow']
            O2flow = H2O_data['O2flow']
            N2flow = H2O_data['N2flow']
            T_cel = H2O_data['T']
            Q1 = H2O_data['Q1']
            Q2 = H2O_data['Q1'] + H2O_data['Q2']
    elif date == 'CLOUD15_NO3_UHL_04Nov22':
            flag_tube = '2'
            R1 = 0.78
            L1 = 10.5
            R2 = 0.5
            L2 = 28.5
            Itx = 5.2009e10 # 6.186e11 * 0.1
            Qx = 20
            file = 'H2O_CLOUD15_04Nov22_NO3_old_attenuator.csv'
            file = input_file_folder + file
            H2O_data = pd.read_csv(file, delimiter=',')
            H2O_1 = H2O_data['H2O_1']
            H2O_2 = H2O_data['H2O_2']
            H2Oconc_1 = H2O_data['H2Oconc_1']
            H2Oconc_2 = H2O_data['H2Oconc_2']
            SO2flow = H2O_data['SO2flow']
            O2flow = H2O_data['O2flow']
            N2flow = H2O_data['N2flow']
            T_cel = H2O_data['T']
            Q1 = H2O_data['Q1']
            Q2 = H2O_data['Q1'] + H2O_data['Q2']
    elif date == 'UFra_comp_test':
            flag_tube = '1'
            R1 = 0.5
            L1 = 39.5
            R2 = 0.5
            L2 = 0
            Itx = 6.186e10
            Qx = 7.6
            file = 'H2O_CLOUD15_UFra_comp_test.csv'
            file = input_file_folder + file
            H2O_data = pd.read_csv(file, delimiter = ',')
            H2O_1 = H2O_data['H2O_1']
            H2O_2 = H2O_data['H2O_2']
            H2Oconc_1 = H2O_data['H2Oconc_1']
            H2Oconc_2 = H2O_data['H2Oconc_2']
            SO2flow = H2O_data['SO2flow']
            O2flow = H2O_data['O2flow']
            N2flow = H2O_data['N2flow']
            T_cel = H2O_data['T']
            Q1 = H2O_data['Q1']
            Q2 = H2O_data['Q1'] + H2O_data['Q2']

    elif date == 'Flotus':
            #%%
            flag_tube = '1'
            R1 = 10 # cm
            L1 = 133.6 # cm
            R2 = 10 # cm
            L2 = 0
            Itx = 5.4e10
            Qx = 20
            file = 'H2O_flotus.csv'
            file = input_file_folder + file
            H2O_data = pd.read_csv(file, delimiter=',')
            H2O_1 = H2O_data['H2O_1']
            H2O_2 = H2O_data['H2O_2']
            H2Oconc_1 = H2O_data['H2Oconc_1']
            H2Oconc_2 = H2O_data['H2Oconc_2']
            SO2flow = H2O_data['SO2flow']
            O2flow = H2O_data['O2flow']
            N2flow = H2O_data['N2flow']
            T_cel = H2O_data['T']
            Q1 = H2O_data['Q1']
            Q2 = H2O_data['Q1'] + H2O_data['Q2']
            #%%
    elif date == 'kinetic_mode':
            flag_tube = '1'
            R1 = 0.5
            L1 = 200
            R2 = 0.5
            L2 = 0
            Itx = 6.186e10
            Qx = 7.6
            file = 'H2O_kinetic_mode.csv'
            file = input_file_folder + file
            H2O_data = pd.read_csv(file, delimiter = ',')
            H2O_1 = H2O_data['H2O_1']
            H2O_2 = H2O_data['H2O_2']
            H2Oconc_1 = H2O_data['H2Oconc_1']
            H2Oconc_2 = H2O_data['H2Oconc_2']
            SO2flow = H2O_data['SO2flow']
            O2flow = H2O_data['O2flow']
            N2flow = H2O_data['N2flow']
            T_cel = H2O_data['T']
            Q1 = H2O_data['Q1']
            Q2 = H2O_data['Q1'] + H2O_data['Q2']
    return R1, L1, R2, L2, flag_tube, file, H2O_1, H2O_2, H2Oconc_1, H2Oconc_2, Q1, Q2, Itx, Qx, SO2flow, O2flow, N2flow, T_cel
#%%
##--------/* Input settings for experiment(s) */--------
## The number of values filled must be the same as that of experiment(s), even values are same
## For example, if two experiments have same temeprature: 'T': [30.5, 30.5]+273.15
para={'P': np.array([101000], dtype=np.float64), #Pressure, Pa
      'outflowLocation': ['before'], # outflow tube located 'before' or 'after' injecting air, water, and so2
      'fullOrSimpleModel': ['full'], # 'simple': Gormley&Kennedy approximation, 'full': flow model (much slower)
      'sampleflow': np.array([30],dtype=np.float64), # inlet sample flow of CIMs, lpm
      'SO2ratio': np.array([1E0],dtype=np.float64)*1e-6, # SO2 ratio of the gas bottle, in ppb
      'O2ratio': np.array([0.209],dtype=np.float64), # O2 ratio in synthetic air
      'Zgrid_num': np.array([80],dtype=np.float64), # number of grids in direction of tube length
      'Rgrid_num': np.array([80],dtype=np.float64), # number of grids in direction of radius
      'dt': np.array([1e-3],dtype=np.float64), # Differential time interval,
      'model_mode': 'normal', #use 'normal' if you don't know what this is for. 'kinetic' mode refers to running
                               #the model without chemistry module to test the kinetic core.
      'date': ['Flotus'], # Experiment date(s) which should be the same as the one(s) in function inputs_setup above
      } 


#%%
##--------/* Run flowtube model */--------
#from Run_flowtube import Run_flowtube
#Run_flowtube(para, folder_flowtube, export_file_folder, inputs_setup) # This is also recorded in the .csv file under export folder

from files_can_deleted.Run_flowtube_for_flotus import Run_flowtube
Run_flowtube(para, folder_flowtube, export_file_folder, inputs_setup)