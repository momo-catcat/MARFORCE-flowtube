# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 15:59:02 2021

@author: jiali
"""

#%% SA_model_mian
import os
import pandas as pd
import numpy as np
from cmd_calib3 import cmd_calib3
from Funcs.Vapour_calc import H2O_conc as H2O_conc


#%% 1. SA_cali model at 2021.9.10
path = os.getcwd()

file=os.path.join(path, "input_files/H2O_1.csv")

H2O_data=pd.read_csv(file)
# delete some rows
H2O_data=H2O_data.drop(index=[0,1,2,3,15])
H2O_data=H2O_data.reset_index()
H2O_data=H2O_data.drop(columns='index')
H2O_data.columns=['order','water','UV']

'''
# setup: box + 41 cm 3/4 inch +58.5 cm 1 inch
# first we roughfuly calculate with 3/4 inch at 99.5 cm long

3/4 inch inner R is 0.78 cm
1 inch inner diamerters is 0.78/4*3
is 2.54 cm, 0.2 for the tube wall
'''

T0 = 273.15
T_cel = 25 #k
T = T0 + T_cel # K
p = 101000 * 1.005 # Pa

ID = 0.78*10*2 # mm the inner diamaters of the tube 
L = 41 * 10 # mm
Q = 22 # lpm

Itx = 5.2009e10 # at Qx flow rate
Qx = 20 # lpm

N2Flow = 22 # slpm
AirFlow = 50 #  synthetic air slpm
O2inAir = 0.209
WaterFlow = H2O_data['water'] #np.array([50, 100, 200]) #

SO2Flow = 5 # slpm
SO2BottlePpm = 5000 # ppm

outflowLocation = 'before' # outflow tube located before or after injecting air, water, and so2

fullOrSimpleModel = 'full' # simple: Gormley&Kennedy approximation, full: flow model (much slower)


#flow time
time = L / (Q * 1e6 / 60 / np.pi / (ID / 2) ** 2)


# computation begins

if outflowLocation in 'after':
    totFlow = N2Flow + AirFlow / 1000 + WaterFlow / 1000 + SO2Flow / 1000
else:
    print(WaterFlow.shape)
    totFlow = Q * np.ones(WaterFlow.shape)

O2conc = O2inAir * AirFlow / 1000 / totFlow * p / 1.3806488e-23 / T / 1e6

H2Oconc = WaterFlow / 1000 / totFlow * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6

SO2conc = SO2Flow / 1000 / totFlow * SO2BottlePpm * 1e-6 * p / 1.3806488e-23 / T / 1e6

It = Itx * Qx / Q

ID1 = 0.78*10*2/3*4
L1 = 68 *10

H2Oconc_1 = H2Oconc
H2SO4 = np.zeros(WaterFlow.shape)
Ho2 = np.zeros(WaterFlow.shape)
for i in range(H2SO4.size):#range(H2SO4.size):
    H2SO4[i],Ho2[i],c = cmd_calib3(O2conc[i], H2Oconc[i], H2Oconc_1[i], SO2conc[i], ID / 10 / 2, L / 10, Q * 1000 / 60,  ID1 / 10 / 2, L1 / 10, Q * 1000 / 60, It, T, p, fullOrSimpleModel, time)


print(H2SO4)

model_1=pd.DataFrame(H2SO4)
model_1['HO2'] =Ho2
model_1.columns= ['SA','HO2']
file=os.path.join(path, "input_files/SA_model_1.csv")
model_1.to_csv(file)
#%% 2. SA_cali model again at 2021.10.28
 
import os
import pandas as pd
import numpy as np
from Funcs.cmd_calib5 import cmd_calib5
path = os.getcwd()

file=os.path.join(path, "input_files/H2O_2.csv")

H2O_data=pd.read_csv(file)
# delete some rows
H2O_data=H2O_data.drop(index=[0,1])
H2O_data=H2O_data.reset_index()
H2O_data=H2O_data.drop(columns='index')
H2O_data.columns=['order','H2O_free','H2O_cali','UV']

'''
# setup: box + 50 cm 3/4 inch +68 cm 1 inch+tower
# first we roughfuly calculate with 3/4 inch at 118 cm long

3/4 inch inner R is 0.78 cm
1 inch inner diamerters is 0.78/4*3
is 2.54 cm, 0.2 for the tube wall
'''
# 50*0.78**2*np.pi/(Q *1000/60)+68*1.04**2*np.pi/(Q1 *1000/60)
T0 = 273.15
T_cel = 25 #k
T = T0 + T_cel # K
p = 101000 * 1.005 # Pa

ID = 0.78*10*2 # mm the inner diamaters of the tube 
L = 50 * 10 # mm
Q = 11 # lpm

Itx = 5.2009e10 # at Qx flow rate
Qx = 20 # lpm

N2Flow = 22 # slpm
AirFlow = 50 #  synthetic air slpm
O2inAir = 0.209
WaterFlow =np.array([100]) #H2O_data['H2O_cali'] #np.array([50, 100, 200]) #[13:19].reset_index()

SO2Flow = 5 # slpm
SO2BottlePpm = 5000 # ppm

outflowLocation = 'before' # outflow tube located before or after injecting air, water, and so2

fullOrSimpleModel = 'full' # simple: Gormley&Kennedy approximation, full: flow model (much slower)


#flow time
time = L / (Q * 1e6 / 60 / np.pi / (ID / 2) ** 2)


# computation begins

if outflowLocation in 'after':
    totFlow = N2Flow + AirFlow / 1000 + WaterFlow / 1000 + SO2Flow / 1000
else:
    print(WaterFlow.shape)
    totFlow = Q * np.ones(WaterFlow.shape)

O2conc = O2inAir * AirFlow / 1000 / totFlow * p / 1.3806488e-23 / T / 1e6

H2Oconc = WaterFlow / 1000 / totFlow * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6

SO2conc = SO2Flow / 1000 / totFlow * SO2BottlePpm * 1e-6 * p / 1.3806488e-23 / T / 1e6

It = Itx * Qx / Q
# second tube 
ID1 = 0.78*10*2/3*4
L1 = 68 *10
Q1 = 23

WaterFlow1 =np.array([0])    #H2O_data['H2O_free']

totFlow1 =Q1 * np.ones(WaterFlow.shape)

H2Oconc_1 =WaterFlow1 / 1000 / totFlow1 * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6

H2SO4 = np.zeros(WaterFlow.shape)
Ho2 = np.zeros(WaterFlow.shape)
for i in range(H2SO4.size):#range(H2SO4.size):
    H2SO4[i],Ho2[i]= cmd_calib5(O2conc[i], H2Oconc[i], H2Oconc_1[i], SO2conc[i], ID / 10 / 2, L / 10, Q * 1000 / 60, ID1 / 10 / 2, L1 / 10, Q1 * 1000 / 60, It, T, p, fullOrSimpleModel, time)
    # H2SO4[i],Ho2[i] = cmd_calib2(O2conc[i], H2Oconc[i], H2Oconc_1[i], SO2conc[i], ID / 10 / 2, L / 10, Q * 1000 / 60, ID1 / 10 / 2, L1 / 10, Q1 * 1000 / 60, It, T, p, fullOrSimpleModel, time )

print(H2SO4)
# h2so4_2 = H2SO4
# h2so4_1 = H2SO4
# model_2=pd.DataFrame(H2SO4)
# model_2['HO2'] =Ho2
# model_2.columns= ['SA','HO2']
# file=os.path.join(path, "input_files/SA_model_2_3.csv") # 2_2 the results from the newest version 09:41-17:00
# model_2.to_csv(file)

#%% 3. SA_cali model again at 2021.11.18
path = os.getcwd()

file=os.path.join(path, "input_files/H2O_3.csv")

H2O_data=pd.read_csv(file)
# delete some rows
H2O_data=H2O_data.drop(index=[0,1,2])
H2O_data=H2O_data.reset_index()
H2O_data=H2O_data.drop(columns='index')
H2O_data.columns=['order','H2O_free','H2O_cali','UV']

'''
# setup: box + 50 cm 3/4 inch +66 cm 1 inch+tower + 21 cm for 1 inch 
# first we roughfuly calculate with 3/4 inch at 116 cm long

3/4 inch inner R is 0.78 cm
1 inch inner diamerters is 0.78/4*3
is 2.54 cm, 0.2 for the tube wall
'''

T0 = 273.15
T_cel = 25 #k
T = T0 + T_cel # K
p = 101000 * 1.005 # Pa

ID = 0.78*10*2 # mm the inner diamaters of the tube 
L = 50 * 10 # mm
Q = 11 # lpm

Itx = 5.2009e10 # at Qx flow rate
Qx = 20 # lpm

N2Flow = 22 # slpm
AirFlow = 50 #  synthetic air slpm
O2inAir = 0.209
WaterFlow = H2O_data['H2O_cali'] #np.array([50, 100, 200]) #

SO2Flow = 5 # slpm
SO2BottlePpm = 5000 # ppm

outflowLocation = 'before' # outflow tube located before or after injecting air, water, and so2

fullOrSimpleModel = 'full' # simple: Gormley&Kennedy approximation, full: flow model (much slower)


#flow time
time = L / (Q * 1e6 / 60 / np.pi / (ID / 2) ** 2)


# computation begins

if outflowLocation in 'after':
    totFlow = N2Flow + AirFlow / 1000 + WaterFlow / 1000 + SO2Flow / 1000
else:
    print(WaterFlow.shape)
    totFlow = Q * np.ones(WaterFlow.shape)

O2conc = O2inAir * AirFlow / 1000 / totFlow * p / 1.3806488e-23 / T / 1e6

H2Oconc = WaterFlow / 1000 / totFlow * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6

SO2conc = SO2Flow / 1000 / totFlow * SO2BottlePpm * 1e-6 * p / 1.3806488e-23 / T / 1e6

It = Itx * Qx / Q

# second tube 
ID1 = 0.78*10*2/3*4
L1 = 66 *10
Q1 = 22
WaterFlow1 = H2O_data['H2O_free']

totFlow1 = Q1 * np.ones(WaterFlow.shape)

H2Oconc_1 = WaterFlow1 / 1000 / totFlow1 * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6

H2SO4 = np.zeros(WaterFlow.shape)
Ho2 = np.zeros(WaterFlow.shape)
for i in range(H2SO4.size):#range(H2SO4.size):
    H2SO4[i],Ho2[i] = cmd_calib3(O2conc[i], H2Oconc[i], H2Oconc_1[i], SO2conc[i], ID / 10 / 2, L / 10, Q * 1000 / 60, ID1 / 10 / 2, L1 / 10, Q1 * 1000 / 60, It, T, p, fullOrSimpleModel, time )

print(H2SO4)
model_3=pd.DataFrame(H2SO4)
model_3['HO2'] =Ho2
model_3.columns= ['SA','HO2']
file=os.path.join(path, "input_files/SA_model_3.csv")
model_3.to_csv(file)