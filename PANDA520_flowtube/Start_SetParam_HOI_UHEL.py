import numpy as np
from Calcu_by_flow_HOI import calculate_concs

'''''''''
set parameters
'''''''''

## 13Sep21. HOI calibration with straight line
## 25Oct21. HOI calibration with Y piece at tower 1

date = '13Sep21'
# date = '25Oct21'
# date = '20Nov21'

if date == '13Sep21':
    paras = dict(p=np.array(101000, dtype=np.float64),  # Pressure, Pa
                 T=np.array(298, dtype=np.float64),  # Temperature, K
                 R1=np.array(0.78, dtype=np.float64),  # ID for the 1st tube
                 L1=np.array(41, dtype=np.float64),  # length for the 1st tube
                 R2=np.array(1.04, dtype=np.float64),  # ID for the 1st tube
                 L2=np.array(58.5, dtype=np.float64),  # length for the 1st tube
                 Itx=5.42e10,  # add it product at the Qx, if you don't have it then you need to calculate it
                 Qx=20,  # Qx where it product was calcuated
                 outflowLocation='before',  # outflow tube located 'before' or 'after' injecting air, water, and so2
                 fullOrSimpleModel='full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (much slower)
                 sampleflow=np.array(22.5, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                 O2ratio=np.array(0.209, dtype=np.float64),  # O2 ratio in synthetic air
                 Zgrid=40,  # number of grids in direction of tube length, usually we use 80
                 Rgrid=80,  # number of grids in direction of radius, usually we use 40
                 dt=np.array(1e-4, dtype=np.float64),  # Differential time interval
                 model_mode='normal',
                 # use 'normal' if you don't know what this is for. 'kinetic' mode refers to running the model without chemistry module to test the kinetic core.
                 Diff_setname=['OH', 'HO2'],
                 # diffusion for the species that you want to define by yourself, otherwise it will be calculated automatically based on the elements it contains
                 Diff_set=[0.215, 0.141],  # add the value according to the Diff_setname
                 sch_name='HOI_cali_chem.txt',  # chemical scheme file name store in the input_mechanism folder
                 const_comp=['I2', 'O2', 'H2O'],
                 # Species you think they should have constant concentrations in the whole tube
                 Init_comp=['OH', 'HO2'],
                 # Species you think they should have initial concentration in the first frid of tube
                 key_spe_for_plot='HOI',  # key species as criterion to stop the loop
                 plot_spec=['OH', 'HOI', 'HO2', 'I', 'I2'],  # species that you want to plot
                 file_name='HOI_cali_T1_13Sep21.csv',
                 # the file you store all the flow data including N2, O2, SO2 in the folder input files
                 )
elif date == '25Oct21':
    paras = dict(p=np.array(101000, dtype=np.float64),  # Pressure, Pa
                 T=np.array(298, dtype=np.float64),  # Temperature, K
                 R1=np.array(0.78, dtype=np.float64),  # ID for the 1st tube
                 L1=np.array(50, dtype=np.float64),  # length for the 1st tube
                 R2=np.array(1.04, dtype=np.float64),  # ID for the 1st tube
                 L2=np.array(68, dtype=np.float64),  # length for the 1st tube
                 Itx=5.42e10,  # add it product at the Qx, if you don't have it then you need to calculate it
                 Qx=20,  # Qx where it product was calcuated
                 outflowLocation='before',  # outflow tube located 'before' or 'after' injecting air, water, and so2
                 fullOrSimpleModel='full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (much slower)
                 sampleflow=np.array(22.5, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                 O2ratio=np.array(0.209, dtype=np.float64),  # O2 ratio in synthetic air
                 Zgrid=40,  # number of grids in direction of tube length, usually we use 80
                 Rgrid=80,  # number of grids in direction of radius, usually we use 40
                 dt=np.array(1e-4, dtype=np.float64),  # Differential time interval
                 model_mode='normal',
                 # use 'normal' if you don't know what this is for. 'kinetic' mode refers to running the model without chemistry module to test the kinetic core.
                 Diff_setname=['OH', 'HO2'],
                 # diffusion for the species that you want to define by yourself, otherwise it will be calculated automatically based on the elements it contains
                 Diff_set=[0.215, 0.141],  # add the value according to the Diff_setname
                 sch_name='HOI_cali_chem.txt',  # chemical scheme file name store in the input_mechanism folder
                 const_comp=['I2', 'O2', 'H2O'],
                 # Species you think they should have constant concentrations in the whole tube
                 Init_comp=['OH', 'HO2'],
                 # Species you think they should have initial concentration in the first frid of tube
                 key_spe_for_plot='HOI',  # key species as criterion to stop the loop

                 plot_spec=['OH', 'HOI', 'HO2', 'I', 'I2'],  # species that you want to plot
                 file_name='HOI_cali_T1_25Oct21.csv',
                 # the file you store all the flow data including N2, O2, SO2 in the folder input files
                 )
elif date == '20Nov21':
    paras = dict(p=np.array(101000, dtype=np.float64),  # Pressure, Pa
                 T=np.array(298, dtype=np.float64),  # Temperature, K
                 R1=np.array(0.78, dtype=np.float64),  # ID for the 1st tube
                 L1=np.array(50, dtype=np.float64),  # length for the 1st tube
                 R2=np.array(1.04, dtype=np.float64),  # ID for the 1st tube
                 L2=np.array(66, dtype=np.float64),  # length for the 1st tube
                 Itx=5.42e10,  # add it product at the Qx, if you don't have it then you need to calculate it
                 Qx=20,  # Qx where it product was calcuated
                 outflowLocation='before',  # outflow tube located 'before' or 'after' injecting air, water, and so2
                 fullOrSimpleModel='full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (much slower)
                 sampleflow=np.array(22.5, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                 O2ratio=np.array(0.209, dtype=np.float64),  # O2 ratio in synthetic air
                 Zgrid=40,  # number of grids in direction of tube length, usually we use 80
                 Rgrid=80,  # number of grids in direction of radius, usually we use 40
                 dt=np.array(1e-4, dtype=np.float64),  # Differential time interval
                 model_mode='normal',
                 # use 'normal' if you don't know what this is for. 'kinetic' mode refers to running the model without chemistry module to test the kinetic core.
                 Diff_setname=['OH', 'HO2'],
                 # diffusion for the species that you want to define by yourself, otherwise it will be calculated automatically based on the elements it contains
                 Diff_set=[0.215, 0.141],  # add the value according to the Diff_setname
                 sch_name='HOI_cali_chem.txt',  # chemical scheme file name store in the input_mechanism folder
                 const_comp=['I2', 'O2', 'H2O'],
                 # Species you think they should have constant concentrations in the whole tube
                 Init_comp=['OH', 'HO2'],
                 # Species you think they should have initial concentration in the first frid of tube
                 key_spe_for_plot='HOI',  # key species as criterion to stop the loop

                 plot_spec=['OH', 'HOI', 'HO2', 'I', 'I2'],  # species that you want to plot
                 file_name='HOI_cali_T2_20Nov21.csv',
                 # the file you store all the flow data including N2, O2, SO2 in the folder input files
                 )


#%%
'''''''''
Calculate the input concentrations based on the parameters
This case only have one tube (no convert) and H2O concentration are unknown 

if you have H2Oconc 
'''''''''

O2conc, I2conc, H2Oconc, paras, export_file_folder = calculate_concs(paras)

'''''''''
prepare the input concentration
'''''''''

if paras['model_mode'] == 'kinetic':  # This mode runs the code in kinetic mode in which chemistry does not exist
    # define the H2SO4 concentration as 1e8 for convenience. !!!! This needs improvement
    OHconc = paras['OHconc'] * 0 + 1e8
    paras['Init_comp'] = ['OH', 'HO2',
                          'H2SO4']  # here you need to change the initial compounds that you already set in paras
    Init_comp_conc = np.transpose([OHconc, OHconc, OHconc])
else:
    # set the initial concentrations for species that you already set in the paras
    Init_comp_conc = np.transpose([paras['OHconc'], paras['OHconc']])

const_comp_conc = np.transpose([I2conc, O2conc, H2Oconc])

## in this case, the number of stages refer to OHconc since we want to check the stage without light on
## in other case, this number need to be added by user
num_stage = paras['OHconc']


# %%
##--------/* Run flowtube model */--------

from Run_flowtube_simplified import Run_flowtube

### to run the flowtube, you need input the const
Run_flowtube(paras, export_file_folder, const_comp_conc, Init_comp_conc, num_stage)
