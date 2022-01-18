# this file is used for all the inputs 
#%% import functions
import numpy as np
import pandas as pd
from Vapour_calc import H2O_conc as H2O_conc
from cmd_calib5 import cmd_calib5
import pene_rate
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
#% Prepare the inputs 

T_cel = 25  # C
T = T_cel + 273.15 # K
p = 101000 * 1.005 # pressure Pa

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

#% set the parameters for the first tube 
 
R1 = 0.78 # cm the inner diameters of the tube
L1 = 50  # cm
Q1 = 22 # lpm

# load H2O 
file= "C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/H2O_1.csv"
H2O_data=pd.read_csv(file)

WaterFlow1 =H2O_data['H20/N2 (mlpm)'][4:-1]

#flow time
# time = L1 / (Q1 * 1e6 / 60 / np.pi / (ID1 / 2) ** 2)

# calculate the conc for SO2, O2, H2O 
if outflowLocation in 'after':
    totFlow = N2Flow + AirFlow / 1000 + WaterFlow1 / 1000 + SO2Flow / 1000
else:
    print(WaterFlow1.shape)
    totFlow = Q1 * np.ones(WaterFlow1.shape)

O2conc1 = O2inAir * AirFlow / 1000 / totFlow * p / 1.3806488e-23 / T / 1e6
H2Oconc1 = WaterFlow1 / 1000 / totFlow * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6
SO2conc1 = SO2Flow / 1000 / totFlow * SO2BottlePpm * 1e-6 * p / 1.3806488e-23 / T / 1e6

It = Itx * Qx / Q1

#% set the parameters for the second tube 
# if there is no second tube, then set to 0

R2 = 0.78/3*4

L2 = 68 

Q2 =  22

WaterFlow2 =np.array([0])    # second H2O flow 

totFlow1 =Q2 * np.ones(WaterFlow1.shape)

H2Oconc2 =WaterFlow2 / 1000 / totFlow1 * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6
O2conc2 =  O2conc1 * Q1/Q2
SO2conc2 = SO2conc1* Q1/Q2

#% save the const comp conc to one file use the same order as const_comp = ['SO2','H2O','O2']
if R1 == 0:
    H2Oconc2 = 0

H2Oconc = np.transpose([H2Oconc1,H2Oconc2])
O2conc =  np.transpose([O2conc1,O2conc2])
SO2conc =  np.transpose([SO2conc1,SO2conc2])

const_comp_conc = np.transpose([SO2conc,H2Oconc,O2conc])

# Initial OH concentration
csH2O = 7.22e-20 #cm2

qyH2O = 1

OHconc = It * csH2O * qyH2O * H2Oconc1
OHconc = OHconc.reset_index()['H20/N2 (mlpm)']

Init_comp_conc = np.transpose([OHconc,OHconc])


rh =  H2Oconc1 * 1e6 * 1.3806488e-23 * T / H2O_conc(T - 273.15, 1).SatP[0]
rh = rh.reset_index()['H20/N2 (mlpm)']


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

# can not find the source use the old one
dHSO3 = 0.126  # cm2/s 
T0 = 300
dHSO3 = 101325 / p * dHSO3 * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

# # https://www.engineeringtoolbox.com/diffusion-coefficients-d_1404.html gas in water seems not good 
# dSO2 = 1.62e-5  # cm2/s 20C 
# T0 = 273.15 + 20
# dSO2 = 101325 / p * dSO2 * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

# dO2 = 2.01e-5  # cm2/s 20C 
# T0 = 273.15 + 20
# dO2 = 101325 / p * dO2 * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

# ['OH', 'SO2', 'HSO3', 'HO2', 'H2O', 'O2', 'H2O2', 'SO3', 'SA']
D = np.array([dOH, pene_rate.cal_diffu('SO2',p,T), dHSO3, dHO2, dH2O, \
              pene_rate.cal_diffu('O2',p,T), pene_rate.cal_diffu('H2O2',p,T), dHSO3,\
                  pene_rate.cal_diffu('H2SO4',p,T)])

# grid for length and radius directions
Zgrid = np.array(40).astype(int)           # number of grid points in tube length direction
Rgrid = np.array(80).astype(int)           # number of grid points in tube radius direction

# chemistry part 
sch_name = os.getcwd()+ '/input_mechanism/SO2_SA.txt'
chm_sch_mrk = ['{', 'RO2', '+', '', '', ';', '+', ';', '$', '{', ':', ';','}']
formula = ['OH','$\mathdefault{HSO_3}$','$\mathdefault{HO_2}$','$\mathdefault{SO_3}$','SA']
key_spe_for_plot = 'H2SO4'
plot_spec = ['OH','HSO3','HO2','SO3','H2SO4'] # plot species 
Init_comp = ['OH','HO2'] # species have inital concentration 
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
          'D': D, # diffusion 
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
          'formula': formula, # the formula for the plots
          'key_spe_for_plot' : key_spe_for_plot, # key species for ploting 
          'plot_spec' : plot_spec # plot species 
          } 
  
    

for i in range(8):
    print(list(params.keys())[i], list(params.values())[i], list(params.values())[i].unit)
    
    # print(list(params.items())[i],list(params.values())[i].unit) 

#%% computation begins
H2SO4 = np.zeros(WaterFlow1.shape)

HO2 = np.zeros(WaterFlow1.shape)

for i in range(H2SO4.size):#range(H2SO4.size):
    meanConc,c= cmd_calib5(const_comp_conc[:,i,:], params, rh[i], Init_comp_conc[i])
  
print(H2SO4)
#%% save the modelled SA, HO2
model = pd.DataFrame([])
model['H2SO4'] = meanConc[8]
model['HO2'] = meanConc[3]

model.to_csv('C:/Users/jiali/MION2-AMT-paper/MION2-AMT-paper/script/SA_cali/input_files/SA_model_1.csv')

# [4.63329389e+06 5.41992198e+07 2.70230841e+07 4.06223710e+07
#  5.41992198e+07 6.77645075e+07 8.13219526e+07 9.48731668e+07
#  1.08418962e+08 8.13219526e+07 2.70230841e+07]
#%% check the specifiy species
import plot_species
plot_species.plot(c,L1,L2,R1,R2,'SA')
