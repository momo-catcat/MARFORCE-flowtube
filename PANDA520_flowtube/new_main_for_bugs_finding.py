# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 20:33:47 2022

@author: jiali
"""
#%% parameters (that can be changed by the user)
#import packages
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

R= 0.78

L= 55

Q=11*1000/60

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


import numpy as np
import os 
import pickle
import def_mod_var
import chem_sch_SMILES
import sch_interr
import eqn_interr
import eqn_pars
# import scipy.constants as si
# import rrc_calc  
# import rate_coeffs 
import water_calc
import init_conc
import RO2_indices  
import write_rate_file
from Vapour_calc import H2O_conc
import matplotlib.pyplot as plt
from scipy import interpolate
from odesolve3 import odesolve as odesolve



# save all the parameters in to the pickle file 
def_mod_var.def_mod_var(0)

# load the pick file 
input_by_sim = os.getcwd()+'/pickle.pkl'

with open(input_by_sim, 'rb') as pk:
   [sav_nam, sch_name, chm_sch_mrk, xml_name, inname, update_stp, 
   		tot_time, comp0, y0, temp, tempt, RH, RHt, Press,  
   		save_step, const_comp, Compt, injectt, Ct,  
   		con_infl_nam, con_infl_t, con_infl_C, 
   		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
   		accom_comp, accom_val, uman_up, int_tol, dil_fac, partit_cutoff,drh_str, erh_str, testf] = pickle.load(pk) 
pk.close()

f_open_eqn = open(sch_name, mode='r') # open the chemical scheme file
	# read the file and store everything into a list
total_list_eqn = f_open_eqn.readlines()
f_open_eqn.close() # close file
const_comp = ['SO2','O2','H2O']
# 	comp_name = [] # list for chemical scheme names of components in the chemical scheme
# 	comp_smil = [] # list for the SMILE strings of components present in the chemical scheme
comp_name, comp_smil, err_mess, H2Oi=chem_sch_SMILES.chem_scheme_SMILES_extr(sch_name, xml_name, chm_sch_mrk)

eqn_list, aqeqn_list, num_eqn, rrc, rrc_name, RO2_names, eqn_list_on=sch_interr.sch_interr(total_list_eqn, chm_sch_mrk)

[rindx, rstoi, pindx, pstoi, reac_coef, nreac, nprod, jac_stoi, jac_den_indx, njac_g, jac_indx, 				
y_arr, y_rind, uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, rstoi_flat, pstoi_flat, 
rr_arr, rr_arr_p,comp_namelist, comp_list, Pybel_objects,comp_num, RO_indx] = eqn_interr.eqn_interr(num_eqn, eqn_list, aqeqn_list, chm_sch_mrk, comp_name, comp_smil)

# comp_namelist, comp_list
# nreac :  number of the reactant 
# nprod : number of production A + B = C nprod is 1 nreac is 2 
# rstoi_flat : record 1D array of stoichiometries per equation A + B = C is 1 1 1, 

[rindx, pindx, rstoi, pstoi, nreac, nprod, jac_stoi, 
		njac, jac_den_indx, jac_indx, y_arr, y_rind,
		uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, 
		rstoi_flat, pstoi_flat, rr_arr, rr_arr_p, 
		comp_num, RO2_indx, RO_indx,
		HOMRO2_indx, comp_list, 
		Pybel_objects, eqn_num, comp_namelist, 
		comp_name, comp_smil, erf, err_mess, con_C_indx] = eqn_pars.extr_mech(sch_name, chm_sch_mrk, xml_name, 
    comp_name, int_tol,  const_comp,
		drh_str, erh_str, dil_fac, sav_nam)
                                                                    
dydt_trak=comp_namelist

                                                                                  
RO2_indices = RO2_indices.RO2_indices(comp_namelist, RO2_names)                                                                                     

                                                                
                                                                    
if H2Oconc == 0:
    meanWeightedH2SO4 = 0
    

# reaction constants (using cm**3 and s)
kSO2pOH = 1.32e-12 * (T / 300) ** -0.7
kOHpHO2 = 4.8e-11 * np.exp(250 / T)
kOHpOH = 6.9e-31 * (T / 300) ** -0.8 * p / 1.3806488e-23 / T / 1e6
kSO3p2H2O = 3.9e-41 * np.exp(6830.6 / T)
kHSO3pO2 = 1.3e-12 * np.exp(-330 / T)

# Initial OH concentration
csH2O = 7.22e-20 #cm2
qyH2O = 1
#It=1.84e10
OHconc = It * csH2O * qyH2O * H2Oconc

# diffusion constants, [cm**2/s]
# order is: HSO3, SO3, HO2, H2SO4, OH

rh = H2Oconc * 1e6 * 1.3806488e-23 * T / H2O_conc(T - 273.15, 1).SatP[0]
# ['OH', 'SO2', 'HSO3', 'HO2', 'H2O', 'O2', 'H2O2', 'SO3', 'SA']
D = np.array([0.215,0.1, 0.126, 0.141, 0.1, 0.126, 0.1,0.126, 0.08])
# ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']
# D = np.array([0.126, 0.126, 0.141, 0.08, 0.215]) #todo need to calculate Diffusion coefficient properly
# D = np.array([0.4, 0.4, 0.4, 0.4, 0.4])
T0 = np.array([300, 300, 298, 298, 298,298,298,298,298])
D = 101325 / p * D * ((T ** (3 / 2)) / (T0 ** (3 / 2)))
# NA = si.Avogadro
H2O, Psat_water, H2O_mw = water_calc.water_calc(T, rh)

dt = 0.00001                        # timestep [s]
numLoop = 500                      # number of times to run to reach the pinhole of the instrument
timesteps = 10000                    # number of timesteps, dt * timesteps * numLoop is time elapsed in the final solution

Zgrid = np.array(40).astype(int)                         # number of grid points in tube length direction
Rgrid = np.array(80).astype(int)                        # number of grid points in tube radius direction

#Change odd number Rgrid to even number grid
if (Rgrid % 2) != 0:
    Rgrid = Rgrid + 1

#% # do not change
# while np.gcd(Rgrid - 1, 10) != 1:  # make sure we don't get a grid point for r = 0 (would give Inf in calculation)
#     Rgrid = int(Rgrid + 1)
# comp_num  = 5
# dr = np.zeros([int(Rgrid),int(Zgrid),comp_num])
# dx = (L) / (Zgrid - 1)
# dr[:,0:int(Zgrid*L/(L+L1)),:] = 2 * R / (Rgrid - 1)
# dr[:,int(Zgrid*L/(L+L1)):,:] =  2 * R1 / (Rgrid - 1)

# H2Otot = np.zeros([int(Rgrid),int(Zgrid)])
# H2Otot[:,0:int(Zgrid*L/(L+L1))] = H2Oconc
# H2Otot[:,int(Zgrid*L/(L+L1)):] =  H2Oconc

# Qtot = np.zeros([int(Rgrid),int(Zgrid),comp_num])
# Qtot[:,0:int(Zgrid*L/(L+L1)),:] = Q
# Qtot[:,int(Zgrid*L/(L+L1)):,:] =  Q1

# SO2tot = np.zeros([int(Rgrid),int(Zgrid)])
# SO2tot[:,0:int(Rgrid*L/(L+L1))] = SO2conc
# # SO2tot[:,int(Zgrid*L/(L+L1)):] =  SO2conc/2
# SO2tot[:,int(Zgrid*L/(L+L1)):] =  SO2conc
# O2tot = np.zeros([int(Rgrid),int(Zgrid)])
# O2tot[:,0:int(Rgrid*L/(L+L1))] = O2conc
# # O2tot[:,int(Zgrid*L/(L+L1)):] =  O2conc/2
# O2tot[:,int(Zgrid*L/(L+L1)):] =  O2conc
# Rtot = np.zeros([int(Rgrid),int(Zgrid),comp_num])
# Rtot[:,0:int(Rgrid*L/(L+L1)),:] = R
# Rtot[:,int(Zgrid*L/(L+L1)):,:] =  R1
#% initial conditions
# order in c is: HSO3, SO3, HO2, H2SO4, OH
c = np.zeros([Rgrid, Zgrid, comp_num])
c[:, 0, 0] = OHconc  # set [OH] at z = 0
c[:, 0, 3] = OHconc  # set [HO2] at z = 0. This equals OH conc
c[:,:,1] = SO2conc

c[:,:,5] = O2conc
c[:,:,4] = H2Oconc

# c[:,0:int(Rgrid*L/(L+L1)),:] = c [:,0:int(Rgrid*L/(L+L1)),:]
## check and output parameters to file
# D_temp = np.tile(D, (Rgrid * Zgrid))
# D = np.transpose(np.reshape(D_temp, (Rgrid, Zgrid, comp_num)), (0, 1, 2))

# r = np.zeros([int(Rgrid), int(Zgrid),  comp_num])
# r = np.abs(np.array([r[i, :, :] + i - r.shape[0] / 2 + 0.5 for i in range(r.shape[0])])) * dr
write_rate_file.write_rate_file(reac_coef,p, rrc, rrc_name, 0)
Comp0 = comp_namelist
# C0 = [1,1,1,1,1,1,1,1,1]
C0 = c[0,0,:] 
# comp_num  = 9
[y,  y_mw, num_comp, M, y_indx_plot, dydt_vst, 
comp_namelist,  erf, err_mess, NOi, HO2i, NO3i]=init_conc.init_conc(comp_num, Comp0, C0, temp, C0[comp_namelist.index('H2O')], Press, Pybel_objects,
testf, dydt_trak, rindx, pindx, num_eqn[0], nreac, nprod, comp_namelist, Compt, comp_namelist, comp_smil, comp_namelist, RO2_indx, HOMRO2_indx, rstoi, pstoi)                                                                                              
y=c[0,0,:]
import RO2_conc
import rate_coeffs
RO2conc = RO2_conc.RO2_conc(RO2_indices,y)
rate_values, erf, err_mess = rate_coeffs.evaluate_rates(RO2conc, y[comp_namelist.index('H2O')], T, 0, M, M*0.7809, y[comp_namelist.index('O2')], 0,  y[comp_namelist.index('HO2')],0,p)
#%import packages

import numpy as np
import matplotlib.pyplot as plt
initc = c
num = 9


dx = L / (Zgrid - 1)
dr = 2 * R / (Rgrid - 1)
r = np.zeros([Rgrid, Zgrid, 9])
r = np.abs(np.array([r[i, :, :] + i - r.shape[0] / 2 + 0.5 for i in range(r.shape[0])])) * dr

#do calculations in matrix
D_temp = np.tile(D, (Rgrid * Zgrid))
D = np.transpose(np.reshape(D_temp, (Rgrid, Zgrid, num)), (0, 1, 2))




#define production and loss terms
term1 = np.zeros([int(Rgrid), int(Zgrid), num])
term2 = np.zeros([int(Rgrid), int(Zgrid), num])
term3 = np.zeros([int(Rgrid), int(Zgrid), num])

# SO2tot = c[:,:,comp_namelist.index('SO2')]
# O2tot = c[:,:,comp_namelist.index('O2')]
# H2Otot = c[:,:,comp_namelist.index('H2O')]

const_comp = ['SO2','O2','H2O']
u = [0,2,3,6,7,8]
# u = [[0,1,2,3,4,5]
yu =2
#%%
for m in range(yu):
    p_a = - 1. / r[1:Rgrid // 2, 1:-1, :] * (initc[1:Rgrid // 2, 1:-1, :] - initc[0:Rgrid // 2 - 1, 1:-1, :]) / dr
    
    p_b = (initc[2:Rgrid // 2 + 1, 1:-1, :] - 2. * initc[1:Rgrid // 2, 1:-1, :] + initc[0:Rgrid // 2 - 1, 1:-1, :]) / (dr ** 2)
    
    p_c = (initc[1:Rgrid // 2, 2:, :] - 2. * initc[1:Rgrid // 2, 1:-1, :] + initc[1:Rgrid // 2, 0:-2,: ]) / (dx * dx)
    
    term1[1:Rgrid // 2, 1:-1, :] = D[1:Rgrid // 2, 1:-1, :] * (p_a + p_b + p_c) #diffusion of gas molecules
    
    # convection; carried by main flow
    # Refs: 1. https://en.wikipedia.org/wiki/Advection
    #       2. Gormley & Kennedy, 1948, Diffusion from a stream flowing through a cylindrical tube
    term2[1:Rgrid // 2, 1:-1, :] = (2. * Q) / (np.pi * R ** 4) * \
                                      (R ** 2 - r[1:Rgrid // 2, 1:-1, :] ** 2) *\
                        (initc[1:Rgrid // 2, 1:-1, :] - initc[1:Rgrid // 2, 0:-2, :]) / dx # carried by main flow
    
    #calculate the last column (measured by the instrument)
    
    p_a_end = - 1. / r[1:Rgrid // 2, -1, :] * (
                initc[1:Rgrid // 2, -1, :] - initc[0:Rgrid // 2 - 1, -1, :]) / dr
    
    p_b_end = (initc[2:Rgrid // 2 + 1, -1, :] - 2. * initc[1:Rgrid // 2, -1, :] + initc[0:Rgrid // 2 - 1, -1, :]) / (dr ** 2)
    
    p_c_end = (initc[1:Rgrid // 2, -1, :] - 2. * initc[1:Rgrid // 2, -2, :] + initc[1:Rgrid // 2, -2,: ]) / (dx * dx)
    
    term1[1:Rgrid // 2, -1, :] = D[1:Rgrid // 2, -1, :] * (p_a_end + p_b_end + p_c_end)
    
    term2[1:Rgrid // 2, -1, :] = (2. * Q) / (np.pi * R ** 4) * \
                                    (R ** 2 - r[1:Rgrid // 2, -1, :] ** 2) *\
                        (initc[1:Rgrid // 2, -1, :] - initc[1:Rgrid // 2, -2, :]) / dx #carried by main flow
    
    # t = ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']
    # term3[1:Rgrid // 2, 1:, comp_namelist.index('HSO3')] = kSO2pOH *  SO2conc * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] - kHSO3pO2 * initc[1:Rgrid // 2, 1:, comp_namelist.index('HSO3')] * O2conc

    # term3[1:Rgrid // 2, 1:, comp_namelist.index('SO3')] = kHSO3pO2 * O2conc * initc[1:Rgrid // 2, 1:, comp_namelist.index('HSO3')] - kSO3p2H2O * H2Oconc * H2Oconc * initc[1:Rgrid // 2, 1:, comp_namelist.index('SO3')]

    # term3[1:Rgrid // 2, 1:, comp_namelist.index('HO2')] = kHSO3pO2 * O2conc * initc[1:-Rgrid // 2, 1:, comp_namelist.index('HSO3')] - kOHpHO2 * initc[1:Rgrid // 2, 1:, comp_namelist.index('HO2')] * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')]

    # term3[1:Rgrid // 2, 1:, comp_namelist.index('SA')] = kSO3p2H2O * H2Oconc * H2Oconc * initc[1:Rgrid // 2, 1:, comp_namelist.index('SO3')]

    # term3[1:Rgrid // 2, 1:,comp_namelist.index('OH')] = -kSO2pOH * SO2conc * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] - kOHpHO2 * \
    #     initc[1:Rgrid // 2, 1:, comp_namelist.index('HO2')] * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] - 2. \
    #         * kOHpOH * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')]
   
    
   
    
   # # NEW VERSION
    # term3 = np.zeros([int(Rgrid), int(Zgrid), num])
    term3[1:Rgrid // 2, 1:,u] = 0

    
    for comp_na in comp_namelist: # get name of this component
        if comp_na not in const_comp: 
            key_name = str(str(comp_na)+ '_comp_indx')  # get index of this component
            compi = dydt_vst[key_name]
            key_name = str(str(comp_na)+'_res')
            dydt_rec = dydt_vst[key_name]
            key_name = str(str(comp_na) + '_reac_sign')
            reac_sign = dydt_vst[key_name]
            # dydt_for = np.zeros(comp_num)
            reac_count = 0
            for i in dydt_rec[0,:]:
                i = int(i) # ensure reaction index is integer - this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float)
                # print(i)
                
                gprate = initc[1:Rgrid // 2, 1:,rindx[i, 0:nreac[i]]] ** rstoi[i, 0:nreac[i]]
                if (len( rstoi[i, 0:nreac[i]]) > 1): 
                    gprate1= gprate[:,:,0] * gprate[:, :,-1] * rate_values[i]
                    
                else:
                    gprate1= gprate[:,:,0]  * rate_values[i]
                # print(rindx[i, 0:nreac[i]])
                # print(rstoi[i, 0:nreac[i]])
                # print(gprate[5,2,:])

                term3[1:Rgrid // 2, 1:, compi] += reac_sign[reac_count]*(gprate1)
                # print(reac_sign[reac_count])
                reac_count += 1
    comp_na = 'OH'
    # key_name = str(str(comp_na)+ '_comp_indx')  # get index of this component
    # compi = dydt_vst[key_name]
    # oh1 = initc[:,:,compi]
    # ohso3 = initc[:,:,2]
    oh1_term3 = term3[:,:,0]
    # hso31= term3[:,:,compi]
    # gp1 = c[1:-1, 1:,0]
    # gp11 = c[1:-1,1:,1]
    # gp2 = rstoi[i, 0:nreac[i]]
    # gp3 = (gp1**gp2[0])
    # gp33 = (gp11**gp2[1])
    # gp4= gp3[:,:]
    # gp5= gp3[:,:]
    # gp6 = gp3*gp33*rate_values[0]
    # term3[1:-1, 1:, compi] += reac_sign[reac_count]*((gp6))
    # NEW VERSION        
    # c = dt * (term1 - term2 + term3) + initc
    c[0:Rgrid // 2,:,u] = dt * (term1[0:Rgrid // 2,:,u] - term2[0:Rgrid // 2,:,u] + term3[0:Rgrid // 2,:,u]) + initc[0:Rgrid // 2,:,u]
    #c[0:Rgrid // 2, :, u] = dt * (term1[0:Rgrid // 2, :, u] - term2[0:Rgrid // 2, :, u])  + initc[0:Rgrid // 2,:, u]
    
    c[Rgrid // 2:, :, u] = np.flipud(c[0:Rgrid // 2,:,u])
    # c[:,:,1] = SO2tot
    # c[:,:,5] = O2tot
    # c[:,:,4] = H2Otot
    initc = c
    oh1 = initc[:,:,0]
    ohso3 = initc[:,:,2]
    ohso3 = initc[:,:,7]
i = 8
plt.pcolor(np.linspace(0, L, Zgrid),np.linspace(-R, R, Rgrid), c[:, : ,i], shading = 'nearest')
plt.ylabel('SA')
i = 0
plt.pcolor(np.linspace(0, L, Zgrid),np.linspace(-R, R, Rgrid), c[:, : ,i], shading = 'nearest')
plt.ylabel('OH')
#%%
for m in range(yu):
    p_a = - 1. / r[1:Rgrid // 2, 1:-1, :] * (initc[1:Rgrid // 2, 1:-1, :] - initc[0:Rgrid // 2 - 1, 1:-1, :]) / dr
    
    p_b = (initc[2:Rgrid // 2 + 1, 1:-1, :] - 2. * initc[1:Rgrid // 2, 1:-1, :] + initc[0:Rgrid // 2 - 1, 1:-1, :]) / (dr ** 2)
    
    p_c = (initc[1:Rgrid // 2, 2:, :] - 2. * initc[1:Rgrid // 2, 1:-1, :] + initc[1:Rgrid // 2, 0:-2,: ]) / (dx * dx)
    
    term1[1:Rgrid // 2, 1:-1, :] = D[1:Rgrid // 2, 1:-1, :] * (p_a + p_b + p_c) #diffusion of gas molecules
    
    # convection; carried by main flow
    # Refs: 1. https://en.wikipedia.org/wiki/Advection
    #       2. Gormley & Kennedy, 1948, Diffusion from a stream flowing through a cylindrical tube
    term2[1:Rgrid // 2, 1:-1, :] = (2. * Q) / (np.pi * R ** 4) * \
                                      (R ** 2 - r[1:Rgrid // 2, 1:-1, :] ** 2) *\
                        (initc[1:Rgrid // 2, 1:-1, :] - initc[1:Rgrid // 2, 0:-2, :]) / dx # carried by main flow
    
    #calculate the last column (measured by the instrument)
    
    p_a_end = - 1. / r[1:Rgrid // 2, -1, :] * (
                initc[1:Rgrid // 2, -1, :] - initc[0:Rgrid // 2 - 1, -1, :]) / dr
    
    p_b_end = (initc[2:Rgrid // 2 + 1, -1, :] - 2. * initc[1:Rgrid // 2, -1, :] + initc[0:Rgrid // 2 - 1, -1, :]) / (dr ** 2)
    
    p_c_end = (initc[1:Rgrid // 2, -1, :] - 2. * initc[1:Rgrid // 2, -2, :] + initc[1:Rgrid // 2, -2,: ]) / (dx * dx)
    
    term1[1:Rgrid // 2, -1, :] = D[1:Rgrid // 2, -1, :] * (p_a_end + p_b_end + p_c_end)
    
    term2[1:Rgrid // 2, -1, :] = (2. * Q) / (np.pi * R ** 4) * \
                                    (R ** 2 - r[1:Rgrid // 2, -1, :] ** 2) *\
                        (initc[1:Rgrid // 2, -1, :] - initc[1:Rgrid // 2, -2, :]) / dx #carried by main flow
    
    # t = ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']
    term3[1:Rgrid // 2, 1:, comp_namelist.index('HSO3')] = kSO2pOH *  SO2conc * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] - kHSO3pO2 * initc[1:Rgrid // 2, 1:, comp_namelist.index('HSO3')] * O2conc

    term3[1:Rgrid // 2, 1:, comp_namelist.index('SO3')] = kHSO3pO2 * O2conc * initc[1:Rgrid // 2, 1:, comp_namelist.index('HSO3')] - kSO3p2H2O * H2Oconc * H2Oconc * initc[1:Rgrid // 2, 1:, comp_namelist.index('SO3')]

    term3[1:Rgrid // 2, 1:, comp_namelist.index('HO2')] = kHSO3pO2 * O2conc * initc[1:-Rgrid // 2, 1:, comp_namelist.index('HSO3')] - kOHpHO2 * initc[1:Rgrid // 2, 1:, comp_namelist.index('HO2')] * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')]

    term3[1:Rgrid // 2, 1:, comp_namelist.index('SA')] = kSO3p2H2O * H2Oconc * H2Oconc * initc[1:Rgrid // 2, 1:, comp_namelist.index('SO3')]

    term3[1:Rgrid // 2, 1:,comp_namelist.index('OH')] = -kSO2pOH * SO2conc * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] - kOHpHO2 * \
        initc[1:Rgrid // 2, 1:, comp_namelist.index('HO2')] * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] - \
             2 * kOHpOH * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')]
    oh2_term3 = term3[:,:,0]
    gprate2_1 =kSO2pOH * SO2conc * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')]
    gprate2_2 = kOHpHO2 * initc[1:Rgrid // 2, 1:, comp_namelist.index('HO2')] * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')]
    gprate2_3 = 2 * kOHpOH * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')] * initc[1:Rgrid // 2, 1:, comp_namelist.index('OH')]
   # # # NEW VERSION
   #  term3 = np.zeros([int(Rgrid), int(Zgrid), num])
    
   
   #  for comp_na in comp_namelist: # get name of this component
   #      if comp_na not in const_comp: 
   #          key_name = str(str(comp_na)+ '_comp_indx')  # get index of this component
   #          compi = dydt_vst[key_name]
   #          key_name = str(str(comp_na)+'_res')
   #          dydt_rec = dydt_vst[key_name]
   #          key_name = str(str(comp_na) + '_reac_sign')
   #          reac_sign = dydt_vst[key_name]
   #          # dydt_for = np.zeros(comp_num)
   #          reac_count = 0
   #          for i in dydt_rec[0,:]:
   #              i = int(i) # ensure reaction index is integer - this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float)
   #              # print(i)
                
   #              gprate = initc[1:Rgrid // 2, 1:,rindx[i, 0:nreac[i]]] ** rstoi[i, 0:nreac[i]]
   #              # print(rindx[i, 0:nreac[i]])
   #              # print(rstoi[i, 0:nreac[i]])
   #              # print(gprate[5,2,:])
   #              gprate1= gprate[:,:,0] * gprate[:, :,-1] * rate_values[i]
   #              term3[1:Rgrid // 2, 1:, compi] += reac_sign[reac_count]*(gprate1)
   #              # print(reac_sign[reac_count])
   #              reac_count += 1
    # comp_na = 'OH'
    # key_name = str(str(comp_na)+ '_comp_indx')  # get index of this component
    # compi = dydt_vst[key_name]
    # oh1 = initc[:,:,compi]
    # ohso3 = initc[:,:,2]
    
    # hso31= term3[:,:,compi]
    # gp1 = c[1:-1, 1:,0]
    # gp11 = c[1:-1,1:,1]
    # gp2 = rstoi[i, 0:nreac[i]]
    # gp3 = (gp1**gp2[0])
    # gp33 = (gp11**gp2[1])
    # gp4= gp3[:,:]
    # gp5= gp3[:,:]
    # gp6 = gp3*gp33*rate_values[0]
    # term3[1:-1, 1:, compi] += reac_sign[reac_count]*((gp6))
    # NEW VERSION        
    # c = dt * (term1 - term2 + term3) + initc
    c[0:Rgrid // 2,:,u] = dt * (term1[0:Rgrid // 2,:,u] - term2[0:Rgrid // 2,:,u] + term3[0:Rgrid // 2,:,u]) + initc[0:Rgrid // 2,:,u]
    #c[0:Rgrid // 2, :, u] = dt * (term1[0:Rgrid // 2, :, u] - term2[0:Rgrid // 2, :, u])  + initc[0:Rgrid // 2,:, u]
    
    c[Rgrid // 2:, :, u] = np.flipud(c[0:Rgrid // 2,:,u])
    # c[:,:,1] = SO2tot
    # c[:,:,5] = O2tot
    # c[:,:,4] = H2Otot
    initc = c
    initc_2 = c
    oh2 = initc_2[:,:,0]
    ohso3_2 = initc_2[:,:,2]
    ohso3_2 = initc_2[:,:,7]
i = 8
plt.pcolor(np.linspace(0, L, Zgrid),np.linspace(-R, R, Rgrid), c[:, : ,i], shading = 'nearest')
plt.ylabel('SA-old')
i = 0
plt.pcolor(np.linspace(0, L, Zgrid),np.linspace(-R, R, Rgrid), c[:, : ,i], shading = 'nearest')
plt.ylabel('OH-old')