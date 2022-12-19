# this file is used for SA calibration that we had for AMT paper

import numpy as np
from Calcu_by_flow_SA import calculate_concs

'''''''''
set parameters
'''''''''
def find_setup_paras(date):
    if date == 'SA_cali_2021-09-10':
        setup_paras = dict(sampleflow=np.array(22.5, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                   R1=np.array(0.78, dtype=np.float64),  # R for the 1st tube
                   L1=np.array(41, dtype=np.float64),  # length for the 1st tube
                   R2=np.array(1.04, dtype=np.float64),  # R for the 2nd tube
                   L2=np.array(58.5, dtype=np.float64),  # length for the 2nd tube
                   file_name='SA_cali_2021-09-10.csv'  # the file you store all the flow data including N2, O2, SO2 in the folder input files
                   )
    elif date == 'SA_cali_2021-10-28':
        setup_paras = dict(sampleflow=np.array(22.5, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                   R1=np.array(0.78, dtype=np.float64),  # R for the 1st tube
                   L1=np.array(50, dtype=np.float64),  # length for the 1st tube
                   R2=np.array(1.04, dtype=np.float64),  # R for the 2nd tube
                   L2=np.array(68, dtype=np.float64),  # length for the 2nd tube
                   file_name='SA_cali_2021-10-28.csv' # the file you store all the flow data including N2, O2, SO2 in the folder input files
                   )
    elif date == 'SA_cali_2021-11-18':
        setup_paras = dict(sampleflow=np.array(22.5, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                   R1=np.array(0.78, dtype=np.float64),  # R for the 1st tube
                   L1=np.array(50, dtype=np.float64),  # length for the 1st tube
                   R2=np.array(1.04, dtype=np.float64),  # R for the 2nd tube
                   L2=np.array(66, dtype=np.float64),  # length for the 2nd tube
                   file_name='SA_cali_2021-11-18.csv' # the file you store all the flow data including N2, O2, SO2 in the folder input files
                   )
    elif date == 'SA_cali_2022-01-04':
        setup_paras = dict(sampleflow=np.array(22.5, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                           R1=np.array(0.78, dtype=np.float64),  # R for the 1st tube
                           L1=np.array(10, dtype=np.float64),  # length for the 1st tube
                           R2=np.array(1.04, dtype=np.float64),  # R for the 2nd tube
                           L2=np.array(78, dtype=np.float64),  # length for the 2nd tube
                           file_name='SA_cali_2022-01-04.csv'
                           # the file you store all the flow data including N2, O2, SO2 in the folder input files
                           )
    elif date == 'SA_cali_2022-01-27':
        setup_paras = dict(sampleflow=np.array(10.6, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                           R1=np.array(0.78, dtype=np.float64),  # R for the 1st tube
                           L1=np.array(26, dtype=np.float64),  # length for the 1st tube
                           file_name='SA_cali_2022-01-27.csv'
                           # the file you store all the flow data including N2, O2, SO2 in the folder input files
                           )
    elif date == 'SA_cali_2022-02-07':
        setup_paras = dict(sampleflow=np.array(17.6, dtype=np.float64),  # inlet sample flow of CIMS, lpm
                   R1=np.array(0.78, dtype=np.float64),  # R for the 1st tube
                   L1=np.array(10, dtype=np.float64),  # length for the 1st tube
                   R2=np.array(1.04, dtype=np.float64),  # R for the 2nd tube
                   L2=np.array(61, dtype=np.float64),  # length for the 2nd tube
                   file_name='SA_cali_2022-02-07.csv' # the file you store all the flow data including N2, O2, SO2 in the folder input files
                   )
    return setup_paras


setup_paras = find_setup_paras('SA_cali_2021-10-28')
#%%
paras = dict(p=np.array(101000, dtype=np.float64),  # Pressure, Pa
             T=np.array(298, dtype=np.float64),  # Temperature, K
             #R1=np.array(0.78, dtype=np.float64),  # ID for the 1st tube
             #L1=np.array(26, dtype=np.float64),  # length for the 1st tube
             Itx=5.42e10,  # add it product at the Qx, if you don't have it then you need to calculate it
             Qx=20,  # Qx where the it product was calculated
             outflowLocation='before',  # outflow tube located 'before' or 'after' injecting air, water, and so2
             fullOrSimpleModel='full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (much slower)
             #sampleflow=np.array(10.6, dtype=np.float64),  # inlet sample flow of CIMS, lpm
             SO2ratio=np.array(5000, dtype=np.float64) * 1e-6,  # SO2 ratio of the gas bottle, in ppb
             O2ratio=np.array(0.209, dtype=np.float64),  # O2 ratio in synthetic air
             Zgrid=40,  # number of grids in direction of tube length, usually we use 80
             Rgrid=80,  # number of grids in direction of radius, usually we use 40
             dt=np.array(1e-4, dtype=np.float64),  # Differential time interval
             model_mode='normal',
             # use 'normal' if you don't know what this is for. 'kinetic' mode refers to running the model without chemistry module to test the kinetic core.
             Diff_setname=['OH', 'HO2', 'SO3', 'H2SO4'],
             # diffusion for the species that you want to define by yourself, otherwise it will be calculated automatically based on the elements it contains
             Diff_set=[0.215, 0.141, 0.126, 0.088],  # add the value according to the Diff_setname
             sch_name='SO2_SA.txt',  # chemical scheme file name store in the input_mechanism folder
             const_comp=['SO2', 'O2', 'H2O'],
             # Species you think they should have constant concentrations in the whole tube
             Init_comp=['OH', 'HO2'],
             # Species you think they should have initial concentration in the first frid of tube
             key_spe_for_plot='H2SO4',  # key species as criterion to stop the loop
             plot_spec=['OH', 'HSO3', 'HO2', 'SO3', 'H2SO4'],  # species that you want to plot
             #file_name='H2O_5.csv',

             )

paras = paras | setup_paras
#%%

'''''''''
Calculate the input concentrations based on the parameters
This case only have one tube (no convert) and H2O concentration are unknown 

if you have H2Oconc 
'''''''''

O2conc, SO2conc, H2Oconc, paras, export_file_folder = calculate_concs(paras)

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

const_comp_conc = np.transpose([SO2conc, O2conc, H2Oconc])

## in this case, the number of stages refer to OHconc since we want to check the stage without light on
## in other case, this number need to be added by user
num_stage = paras['OHconc']


# %%
##--------/* Run flowtube model */--------

from Run_flowtube_simplified import Run_flowtube

### to run the flowtube, you need input the const
Run_flowtube(paras, export_file_folder, const_comp_conc, Init_comp_conc, num_stage)
