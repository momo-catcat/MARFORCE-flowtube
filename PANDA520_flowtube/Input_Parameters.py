#input packages
#just try to change something in this file
import numpy as np
# from cmd_calib1 import cmd_calib1
from Vapour_calc import H2O_conc as H2O_conc
from cmd_calib5 import cmd_calib5


#calculate vapor pressure


#calculate saturation vapour pressure


# if there are two different tube then go to the cmd_cali1Matlab.m

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

ID = 0.78*10*2 # mm the inner diameters of the tube
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
# ID1 = 0.78*10*2/3*4
ID1 = 0.78*10*2
L1 = 68 *10
# L1 = 0 *10
Q1 =  11

WaterFlow1 =np.array([0])    #H2O_data['H2O_free']

totFlow1 =Q1 * np.ones(WaterFlow.shape)

H2Oconc_1 =WaterFlow1 / 1000 / totFlow1 * H2O_conc(T_cel, 1).SatP[0] / 1.3806488e-23 / T / 1e6

H2SO4 = np.zeros(WaterFlow.shape)
HO2 = np.zeros(WaterFlow.shape)
for i in range(H2SO4.size):#range(H2SO4.size):
    H2SO4[i],HO2[i]= cmd_calib5(O2conc[i], H2Oconc[i], H2Oconc_1[i], SO2conc[i], ID / 10 / 2, L / 10, Q * 1000 / 60, ID1 / 10 / 2, L1 / 10, Q1 * 1000 / 60, It, T, p, fullOrSimpleModel, time)
    # H2SO4[i],Ho2[i] = cmd_calib2(O2conc[i], H2Oconc[i], H2Oconc_1[i], SO2conc[i], ID / 10 / 2, L / 10, Q * 1000 / 60, ID1 / 10 / 2, L1 / 10, Q1 * 1000 / 60, It, T, p, fullOrSimpleModel, time )

print(H2SO4)