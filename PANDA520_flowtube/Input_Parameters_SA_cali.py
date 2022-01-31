# this file is used for all the inputs 
#%% first clear all the variable that you
# already have 
for name in dir():
    if not name.startswith('_'):
        del globals()[name]
del name

#% import functions 
import numpy as np
import pandas as pd
from Vapour_calc import H2O_conc as H2O_conc
from cmd_calib5 import cmd_calib5
import os

class UnitFloat(float):

    def __new__(self, value, unit=None):
       return float.__new__(self, value)

    def __init__(self, value, unit=None):
        self.unit = unit

'''
# setup: box + 50 cm 3/4 inch +68 cm 1 inch+tower
# first we roughfuly calculate with 3/4 inch at 118 cm long
3/4 inch inner R is 0.78 cm
1 inch inner diamerters is 0.78/4*3
is 2.54 cm, 0.2 for the tube wall
'''
#%% Prepare the inputs
  
# load H2O Q  
file = os.getcwd() + '/input_files/test.csv'

H2O_data = pd.read_csv(file)

''' set temperature and press '''
# T_cel = 25  # C

T_cel = float(input('temperaure C:'))
T = T_cel + 273.15 # K
p = 96060 * 1.005 # pressure Pa

#% set the parameters for the first tube 
flag_tube = input('how much tube you have:')

if flag_tube == '2':
    R1 = float(input('1st tube diameter:'))
    L1 = float(input('1st tube length:'))
    R2 = float(input('2nd tube diameter:')) 
    L2 = float(input('2nd tube length:'))
    Q1 = H2O_data['Q1'][0] # lpm
    Q2 =  H2O_data['Q2'][2]
    Q2 = Q1 + Q2
else:
    R1 = float(input('tube diameter:'))
    L1 = float(input('tube length:'))
    R2 = 0 
    L2 = 0
    Q1 = H2O_data['Q1'][0] # lpm
    Q2 = 0
#% set the parameters for the first tube 
# R1 = 0.78 # cm the inner diameters of the tube
# L1 = 10  # cm
# Q1 = H2O_data['Q1'][0] # lpm

#% set the parameters for the second tube 
# if there is no second tube, then set to 0
# R2 = 0.78/3*4
# L2 = 68 
# Q2 =  H2O_data['Q2'][2]
# Q2 = Q1 + Q2

# It product for cali box 
Itx = 5.2009e10 # at Qx flow rate
Qx = 20 # lpm

N2Flow = 22 # slpm
AirFlow = 50 #  synthetic air slpm
O2inAir = 0.209
SO2Flow = 5 # slpm
SO2BottlePpm = 5000 # ppm

outflowLocation = 'before' # outflow tube located before or after injecting air, water, and so2

fullOrSimpleModel = 'full' # simple: Gormley&Kennedy approximation, full: flow model (much slower)

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

WaterFlow2 = H2O_data['H2O_2']   # second H2O flow 
totFlow2 = Q2 * np.ones(WaterFlow1.shape)
H2Oconc2 = (WaterFlow2 + WaterFlow1) / 1000 / totFlow2 * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6
O2conc2 =  O2conc1 * Q1/Q2
SO2conc2 = SO2conc1* Q1/Q2

if flag_tube == '1':
    H2Oconc2 = pd.Series(np.zeros(WaterFlow1.shape))
    O2conc2 = pd.Series(np.zeros(WaterFlow1.shape))
    SO2conc2 = pd.Series(np.zeros(WaterFlow1.shape))

H2Oconc = np.transpose([H2Oconc1,H2Oconc2])
O2conc =  np.transpose([O2conc1,O2conc2])
SO2conc =  np.transpose([SO2conc1,SO2conc2])

#% store all the const species to const_comp_conc
const_comp_conc = np.transpose([SO2conc,H2Oconc,O2conc])

# calculate Initial concentration for some species here is OH and HO2 
csH2O = 7.22e-20 #cm2

qyH2O = 1

It = Itx * Qx / Q1

OHconc = It * csH2O * qyH2O * H2Oconc1
OHconc = OHconc.reset_index()['H2O_1']

# store initial concentration
Init_comp = ['OH','HO2'] # species have inital concentration 
Init_comp_conc = np.transpose([OHconc,OHconc])

# diffusion constants, [cm**2/s]

# DOH-air = 165 ± 20 Torr cm2 s-1, DHO2-He = 430 ± 30 Torr cm2 s-1，DO3-He = 410 ± 25 Torr cm2 s-1 at 296 K. 
#Source OH, HO2, and Ozone Gaseous Diffusion Coefficients
dOH = 165/(0.00750062*p) # convert to pa
T0 = 296
dOH = 101325 / p * dOH * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

#Source A Measurement of the Diffusion Coefficient of Hydrogen Peroxide Vapor into Air
dHO2 = 111/(0.00750062*p) # convert to pa
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

Diff_setname = ['OH','H2O','HO2']
Diff_set = [dOH,dHO2,dH2O]

# grid for length and radius directions
Zgrid = np.array(40).astype(int)           # number of grid points in tube length direction
Rgrid = np.array(80).astype(int)           # number of grid points in tube radius direction

# chemistry part 
sch_name = os.getcwd()+ '/input_mechanism/SO2_SA.txt'
chm_sch_mrk = ['{', 'RO2', '+', '', '', ';', '+', ';', '$', '{', ':', ';','}']
# formula = ['OH','$\mathdefault{HSO_3}$','$\mathdefault{HO_2}$','$\mathdefault{SO_3}$','SA']
key_spe_for_plot = 'H2SO4'
plot_spec = ['OH','HSO3','HO2','SO3','H2SO4'] # plot species 
const_comp = ['SO2','H2O','O2']# species have constant concentration
drh_str = str('0.*TEMP')
erh_str = str('0.*TEMP')

params = {'T' : UnitFloat(T, "K"), # temperaure 
          'p': UnitFloat(p, "Pa"), # pressure pa
          'R1': UnitFloat(R1, "cm"), # diameters for frist tube 
          'R2': UnitFloat(R2, "cm"), # diameters for second tube 
          'L1': UnitFloat(L1, "cm"), # length for frist tube 
          'L2': UnitFloat(L2, "cm"), # length for frist tube 
          'Q1': UnitFloat(Q1, "lpm"), # flow for frist tube 
          'Q2': UnitFloat(Q2, "lpm"), # flow for second tube
          'dt': 0.00001, # flow for second tube
          'Diff_setname': Diff_setname, # diffusion
          'Diff_set': Diff_set,
          # 'It': It, # it product for calculation 
          'fullOrSimpleModel': fullOrSimpleModel, #Gormley&Kennedy approximation, full: flow model (much slower)
          'sch_name': sch_name, # file for the MCM file
          'chm_sch_mrk': chm_sch_mrk,# markers to isolate sections of chemical scheme based on MCM KPP format
          'drh_str':str('0.*TEMP'),
          'erh_str':str('0.*TEMP'),
          'const_comp' : const_comp, # constant species
          'Init_comp': Init_comp, # species have inital concentration 
          'Zgrid': Zgrid, # number of grid points in tube length direction
          'Rgrid': Rgrid, # number of grid points in tube radius direction
          # 'formula': formula, # the formula for the plots
          'key_spe_for_plot' : key_spe_for_plot, # key species for ploting 
          'plot_spec' : plot_spec # plot species 
          } 
#% 
for i in range(8):
    print(list(params.keys())[i], list(params.values())[i], list(params.values())[i].unit)

# i = 6
# const_comp_conc= const_comp_conc[:,i,:]
# Init_comp_conc=Init_comp_conc[i]

#% computation begins
meanconc = []

c = []

for i in range(WaterFlow1.size):#range(H2SO4.size):
    if H2Oconc1[i]>0:  
        meanConc1,c1= cmd_calib5(const_comp_conc[:,i,:], params, Init_comp_conc[i])
        meanconc.append(meanConc1) 
        c.append(c1)

# meanconc = meanConc1
meanconc_s = pd.DataFrame(np.transpose(meanconc)) 
meanconc_s.index = plot_spec
    
#%% save the modelled SA, HO2

meanconc_s = pd.DataFrame(np.transpose(meanconc)) 
meanconc_s.index = plot_spec
meanconc_s.to_csv('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_2_mean.csv')

with open('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_2_c.txt', 'w') as f:
    for item in c:
        f.write("%s\n" % item)



#%% check the specifiy species
import plot_species
plot_species.plot(c,L1,L2,ID1,ID2,'SA') # need to change
