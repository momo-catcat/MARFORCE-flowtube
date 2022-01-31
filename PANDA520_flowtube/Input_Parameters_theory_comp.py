# this file is used for all the inputs 
#%% first clear all the variable that you
# already have 
for name in dir():
    if not name.startswith('_'):
        del globals()[name]
del name

#%% import functions 
import numpy as np
import pandas as pd
from Vapour_calc import H2O_conc as H2O_conc
from cmd_calib_theory_comp import cmd_calib_theory_comp
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
#%%Prepare the inputs
  
''' set temperature and press '''
# T_cel = 25  # C

T_cel = 25
T = T_cel + 273.15 # K
p = 101000 # pressure Pa

#% set the parameters for the first tube 
flag_tube = 2

if flag_tube == 2:
    R1 = 0.78
    L1 = 100
    R2 = 0.78
    L2 = 100
    Q1 = 10
    Q2 = Q1
else:
    R1 = 0.78
    L1 = 200
    R2 = 0 
    L2 = 0
    Q1 = 22 # lpm
    Q2 = 0

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

WaterFlow1 = np.array([50])

if outflowLocation in 'after':
    totFlow = N2Flow + AirFlow / 1000 + WaterFlow1 / 1000 + SO2Flow / 1000
else:
    totFlow = Q1 * np.ones(WaterFlow1.shape)

O2conc1 = O2inAir * AirFlow / 1000 / totFlow * p / 1.3806488e-23 / T / 1e6
H2Oconc1 = WaterFlow1 / 1000 / totFlow * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6
SO2conc1 = SO2Flow / 1000 / totFlow * SO2BottlePpm * 1e-6 * p / 1.3806488e-23 / T / 1e6


H2Oconc = H2Oconc1
O2conc = O2conc1
SO2conc = SO2conc1

#% store all the const species to const_comp_conc
const_comp_conc = np.transpose([SO2conc,H2Oconc,O2conc])
# calculate Initial concentration for some species here is OH and HO2
csH2O = 7.22e-20 #cm2

qyH2O = 1

It = Itx * Qx / Q1

OHconc = 1e8

# store initial concentration
Init_comp = ['H2SO4','OH'] # species have inital concentration
Init_comp_conc = np.transpose([OHconc,OHconc])

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

# https://www.engineeringtoolbox.com/air-diffusion-coefficient-gas-mixture-temperature-d_2010.html
dH2SO4 = 0.0921  # cm2/s 20C

# # can not find the source use the old one
# dHSO3 = 0.126  # cm2/s
# T0 = 300
# dHSO3 = 101325 / p * dHSO3 * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

Diff_setname = ['OH','H2O','HO2','H2SO4']
Diff_set = [dOH,dHO2,dH2O,dH2SO4]

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
for i in range(8):
    print(list(params.keys())[i], list(params.values())[i], list(params.values())[i].unit)
# i =1 
# const_comp_conc= const_comp_conc[:,i,:]
# Init_comp_conc=Init_comp_conc[i]
meanconc = []

c = []

for i in range(WaterFlow1.size):#range(H2SO4.size):
    if H2Oconc1[i]>0:
        meanConc1,c1= cmd_calib_theory_comp(const_comp_conc[0], params, Init_comp_conc[i])
        meanconc.append(meanConc1) 
        c.append(c1)

meanconc_s = pd.DataFrame(np.transpose(meanconc)) 
meanconc_s.index = plot_spec
    



