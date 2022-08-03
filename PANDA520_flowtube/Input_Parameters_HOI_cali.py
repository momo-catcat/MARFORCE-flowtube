# this file is used for all the inputs

# clear variables
for name in dir():
    if not name.startswith('_'):
        del globals()[name]
del name

# %% import packages and functions
import numpy as np
import pandas as pd
import os
from cmd_calib5 import cmd_calib5
from Calcu_by_flow import const_comp_conc_cal, const_comp_conc_cal_H2O, const_comp_conc_cal_OH
from diffusion_const_added import add_diff_const as add_diff_const

# add unit after values
class UnitFloat(float):

    def __new__(self, value, unit=None):
        return float.__new__(self, value)

    def __init__(self, value, unit=None):
        self.unit = unit


# set temperature and pressure
T_cel = 23
T = T_cel + 273.15  # K
p = 101000  # pressure Pa

# load input file for flows and concentrations
file = os.getcwd() + '/input_files/HOI_cali_T1_25Oct21.csv'
#file = os.getcwd() + '/input_files/HOI_cali_T2_20Nov21.csv'
H2O_data = pd.read_csv(file)

''' The format of input file needs to be changed:
H2O_1 and H2O_2 refer to the flows in the 1st and 2nd tube,
I2conc_1 and I2conc_2 refer to the concentrations in the 1st and 2nd tube.
'''

# set flow and input concentrations
Q1 = H2O_data['Q1']
Q2 = H2O_data['Q2']
H2O_1 = H2O_data['H2O_1']
H2O_2 = H2O_data['H2O_2']
I2conc_1 = H2O_data['I2conc_1']
I2conc_2 = H2O_data['I2conc_2']
O2flow = H2O_data['SynAir_cali']

flag_tube = '3'
outflowLocation = 'before'  # outflow tube located before or after injecting air, water, and so2
fullOrSimpleModel = 'full'  # simple: Gormley&Kennedy approximation, full: flow model (much slower)
# in this study, we have the outflow before injecting air, water and I2

# % set the parameters for the first tube
if flag_tube in ['3', '4']:
    # 20Nov21
    L1 = 50
    L2 = 66
    R1 = 0.78
    R2 = 1.04
    Q2 = Q1 + Q2

elif flag_tube == '2':
    R1 = 0.78
    L1 = 100
    R2 = 1.04
    L2 = 50
    Q2 = Q1

else:
    R1 = 0.78
    R2 = 0
    L2 = 0
    L1 = 100
    Q2 = 0

sampflow = 22.5  # lpm
O2ratio = 0.209  # O2inAir = 0.209
N2Flow = 23  # lpm
const_comp_pre = ['H2O', 'O2']  # species have constant concentration and are calculated from flows
const_comp_pre_know = ['I2']  # species have known constant concentration but already known
const_comp = const_comp_pre + const_comp_pre_know  # species have constant concentration
# get all the concentrations
O2conc = const_comp_conc_cal(O2flow, outflowLocation, sampflow, H2O_1, N2Flow, O2ratio,
                             Q1, Q2, T_cel, T, p, flag_tube)

H2Oconc = const_comp_conc_cal_H2O(O2flow, outflowLocation, sampflow, H2O_1, H2O_2, N2Flow, O2ratio,
                                  Q1, Q2, T_cel, T, p, flag_tube)

I2conc = np.transpose([I2conc_1, I2conc_2])

# % store all the const species to const_comp_conc follow the order of const_comp
const_comp_conc = np.transpose([H2Oconc, O2conc, I2conc])

OHconc, const_comp_free, const_comp_conc_free = const_comp_conc_cal_OH(H2Oconc, O2conc, Q1, flag_tube)

# store initial concentration
Init_comp = ['OH', 'HO2']  # species have inital concentration
Init_comp_conc = np.transpose([OHconc, OHconc])

# add more diffusion constants in diffusion_const_added.py file if you want
Diff_setname = ['OH', 'H2O', 'HO2']  # diffusion constants, [cm**2/s]
dOH, dH2O, dHO2 = add_diff_const(p, T)
Diff_set = [dOH, dH2O, dHO2]

# grid for length and radius directions
Zgrid = np.array(40).astype(int)  # number of grid points in tube length direction
Rgrid = np.array(80).astype(int)  # number of grid points in tube radius direction

# chemistry part
sch_name = os.getcwd() + '/input_mechanism/HOI_cali_chem.txt'
chm_sch_mrk = ['{', 'RO2', '+', '', '', ';', '+', ';', '$', '{', ':', ';', '}']
key_spe_for_plot = 'HOI'
plot_spec = ['OH', 'HOI', 'HO2', 'I', 'I2']  # plot species

params = {'T': UnitFloat(T, "K"),  # temperature
          'p': UnitFloat(p, "Pa"),  # pressure pa
          'R1': UnitFloat(R1, "cm"),  # diameters for first tube
          'R2': UnitFloat(R2, "cm"),  # diameters for second tube
          'L1': UnitFloat(L1, "cm"),  # length for first tube
          'L2': UnitFloat(L2, "cm"),  # length for first tube
          'dt': 0.0001,  # flow for second tube
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

for i in range(6):
    print(list(params.keys())[i], list(params.values())[i], list(params.values())[i].unit)

# %% computation begins
meanconc = []

c = []

for i in range(len(H2O_data)):  # range(H2SO4.size):
    if H2O_data['H2O_1'][i] > 0:
        meanConc1, c1 = cmd_calib5(const_comp_conc[:, i, :], params, Init_comp_conc[i], Q1[i], Q2[i])
        meanconc.append(meanConc1)
        c.append(c1)

meanconc_s = pd.DataFrame(meanconc)
meanconc_s.columns = plot_spec
# meanconc_s.to_csv('./Export_files/HOI_cali_25Oct21.csv')
# meanconc_s.to_csv('./Export_files/HOI_cali_20Nov21.csv')

# with open('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_4_c.txt', 'w') as f:
#     for item in c:
#         f.write("%s\n" % item)
