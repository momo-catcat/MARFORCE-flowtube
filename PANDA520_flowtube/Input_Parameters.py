# this file is used for all the inputs 
#%% import functions
import numpy as np
from Vapour_calc import H2O_conc as H2O_conc
from cmd_calib5 import cmd_calib5
import def_mod_var

'''
# setup: box + 50 cm 3/4 inch +68 cm 1 inch+tower
# first we roughfuly calculate with 3/4 inch at 118 cm long

3/4 inch inner R is 0.78 cm
1 inch inner diamerters is 0.78/4*3
is 2.54 cm, 0.2 for the tube wall
'''

#% Prepare the inputs 

T0 = 273.15 # initial 
T_cel = 25 # C
T = T0 + T_cel # K
p = 101000 * 1.005 # Pa

# save all the parameters in to the pickle file change the const_comp
# here we have 'SO2','H2O','O2'

const_comp = ['SO2','H2O','O2']
def_mod_var.def_mod_var(0,T,p,const_comp)  

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
 
ID1 = 0.78 # mm the inner diameters of the tube
L1 = 50  # mm
Q1 = 11 # lpm
WaterFlow1 =np.array([100,500,600,700]) #H2O_data['H2O_cali']

#flow time
time = L1 / (Q1 * 1e6 / 60 / np.pi / (ID1 / 2) ** 2)

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

ID2 = 0.78/3*4

L2 = 68 

Q2 =  22

WaterFlow2 =np.array([100,200,300,500])    # second H2O flow 

totFlow1 =Q2 * np.ones(WaterFlow1.shape)

H2Oconc2 =WaterFlow2 / 1000 / totFlow1 * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6
O2conc2 =  O2conc1 * Q1/Q2
SO2conc2 = SO2conc1* Q1/Q2

#% save the const comp conc to one file use the same order as const_comp = ['SO2','H2O','O2']
if ID1 == 0:
    H2Oconc2 = 0

H2Oconc = np.transpose([H2Oconc1,H2Oconc2])
O2conc =  np.transpose([O2conc1,O2conc2])
SO2conc =  np.transpose([SO2conc1,SO2conc2])

const_comp_conc = np.transpose([SO2conc,H2Oconc,O2conc])

#% computation begins

H2SO4 = np.zeros(WaterFlow1.shape)
HO2 = np.zeros(WaterFlow1.shape)


for i in range(H2SO4.size):#range(H2SO4.size):
    H2SO4[i],HO2[i],c= cmd_calib5(const_comp_conc[:,i,:], ID1, L1, Q1 * 1000 / 60, ID2, L2, Q2 * 1000 / 60, It, fullOrSimpleModel, time)
  
print(H2SO4)
#%% check the specifiy species
import plot_species
plot_species.plot(c,L1,L2,ID1,ID2,'SA')
