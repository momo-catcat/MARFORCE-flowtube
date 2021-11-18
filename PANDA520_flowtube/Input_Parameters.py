#input packages
import numpy as np
from cmd_calib1 import cmd_calib1
from Vapour_calc import H2O_conc as H2O_conc


#calculate vapor pressure


#calculate saturation vapour pressure


# if there are two different tube then go to the cmd_cali1Matlab.m

T0 = 273.15
T_cel = 25 #k
T = T0 + T_cel # K
p = 101000 * 1.005 # Pa

ID = 24 # mm the inner diamaters of the tube
L = 93 * 10 # mm
Q = 22 # lpm

Itx = 5.2009e10 # at Qx flow rate
Qx = 20 # lpm

N2Flow = 22 # slpm
AirFlow = 50 #  synthetic air slpm
O2inAir = 0.209
WaterFlow = np.array([50, 100, 200])

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

H2SO4 = np.zeros(WaterFlow.shape)
for i in range(1):#range(H2SO4.size):
    H2SO4[i] = cmd_calib1(O2conc[i], H2Oconc[i], SO2conc[i], ID / 10 / 2, L / 10, Q * 1000 / 60, It, T, p, fullOrSimpleModel, time)


