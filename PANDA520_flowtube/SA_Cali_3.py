# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 17:38:19 2021

@author: jiali
"""

#%% analysis the data of api_data_3 
#directly load the csv of SA_cali_3 file
import os
import glob
import sys
import pandas as pd
import numpy as np 
import datetime
import scipy.io
import scipy.stats
import types
import matplotlib.pyplot as plt
from read_tofware_preprocessed_files import quick_plot
from detect_limit import get_detect_llimit
path = os.getcwd()

file=os.path.join(path, "input_files/api_data_3.csv")

api_data=pd.read_csv(file)
api_data['time']=pd.to_datetime(api_data['time'])
api_data['time']=api_data['time']- datetime.timedelta(hours=1)
api_data=api_data.set_index('time')

Br=api_data['Br-']

H2O_Br=api_data['H2OBr-']

Primary_ions=Br+H2O_Br

Sa=(api_data['H2SO4Br-'])/Primary_ions

H2O=H2O_Br/Primary_ions

HO2=api_data['HO2Br-']/Primary_ions

SO2=api_data['SO2Br-']/Primary_ions



#% load time table for SA cali

'''
the time table for the experiment and the time from tofware
there is one hour different, the tofware time is more one hour later
'''

file_name= os.path.join(path, "input_files/SA_cali_3.xlsx")

time_table=pd.read_excel(file_name, index_col=None, header=2,engine='openpyxl')

time_table['Start time']=pd.to_datetime(time_table['Start time'])
time_table['End time']=pd.to_datetime(time_table['End time'])
# check the table remove nan rows
print(time_table)

time_table=time_table.drop(index=[24,25,26,27]) 

print(time_table.columns)

#% delete the useless line like test and before 
# time_table=time_table.drop(index=[0,1,4])
# time_table=time_table.reset_index(drop=True)

# check the H2O and UV lights 
H2o=time_table[['H2O free (MFC8)','H2O cali (MFC4+MFC2)','UVC.1 filter']]
# save the water flow to the folder
H2o.to_csv(path+"/input_files/H2O_3.csv")

Sa_stage=[]
H2O_stage=[]
Ho2_stage=[]
for i in range(len(time_table)):
    Sa_stage.append(Sa[time_table['Start time'][i]:time_table['End time'][i]])
    Ho2_stage.append(HO2[time_table['Start time'][i]:time_table['End time'][i]])
    H2O_stage.append(H2O[time_table['Start time'][i]:time_table['End time'][i]])

Sa_mean_API=[Sa_stage[i][3:-3].mean() for i in range(len(Sa_stage))] 
H2O_mean_API=[H2O_stage[i][3:-3].mean() for i in range(len(H2O_stage))] 
Ho2_mean_API=[Ho2_stage[i][3:-3].mean() for i in range(len(Ho2_stage))] 
det_limit_sa=get_detect_llimit(Sa_stage[0])
#% check the SA measured mean values
quick_plot(Sa)
sa=pd.concat([Sa_stage[19],Sa_stage[20]],axis=0)
Sa_mean_API[6]
sa=Sa_stage[22]
len(sa)
sa=Sa_stage[22][30:]  
quick_plot(sa)
'''
stage 6 is 3:-5 
stage 9 is 22:-2 
after stage 10-12 is not good can not find the good data 
stage 10 is 21: 
stage 11 is 20:-1 
stage 12 is 25:-1 
stage 13 is 15:-3
stage 14 is 11:-3
stage 15 is 12:-3
stage 19 is 13:-3
stage 20 is 15:-3
stage 21 is 24:
stage 22 is 38:-2
stage 23 is 30:   
'''
Sa_mean_API[6]=Sa_stage[6][3:-5].mean()
Sa_mean_API[9]=Sa_stage[9][22:-2].mean()
Sa_mean_API[10]=Sa_stage[10][21].mean()
Sa_mean_API[11]=Sa_stage[11][20:-1].mean()
Sa_mean_API[12]=Sa_stage[12][25:-1].mean()
Sa_mean_API[13]=Sa_stage[13][15:-3].mean()
Sa_mean_API[14]=Sa_stage[14][11:-3].mean()
Sa_mean_API[15]=Sa_stage[15][12:-3].mean()
Sa_mean_API[19]=Sa_stage[19][13:-3].mean()
Sa_mean_API[20]=Sa_stage[20][15:-3].mean()
Sa_mean_API[21]=Sa_stage[21][24:].mean()
Sa_mean_API[22]=Sa_stage[22][39:-2].mean()
Sa_mean_API[23]=Sa_stage[23][30:].mean()
Sa_mean_API_3=Sa_mean_API
#% check the H2O
quick_plot(H2O)
h2o=pd.concat([H2O_stage[19],H2O_stage[20]],axis=0)
H2O_mean_API[6]
h2o=H2O_stage[23]
len(h2o)
h2o=H2O_stage[21][25:]  
quick_plot(h2o)

'''
stage 0 is 50: 
stage 6 is 0:-10
stage 9 is 25:
stage 10 is 18:-1
stage 14 is 10:-3
stage 20 is 10:-4
stage 21 is 25:
'''
H2O_mean_API[0]=H2O_stage[0][50:].mean()
H2O_mean_API[6]=H2O_stage[6][0:-10].mean()
H2O_mean_API[9]=H2O_stage[9][25:].mean()
H2O_mean_API[10]=H2O_stage[10][18:-1].mean()
H2O_mean_API[14]=H2O_stage[14][10:-3].mean()
H2O_mean_API[20]=H2O_stage[20][10:-4].mean()
H2O_mean_API[21]=H2O_stage[21][25:].mean()
H2O_stage_conc_3 =[H2O_mean_API[i] *10 ** 18.02663192 for i in range(len(H2O_mean_API))]

#% check the HO2 measured mean values
quick_plot(HO2)
ho2=pd.concat([Ho2_stage[0],Ho2_stage[1]],axis=0)
Ho2_mean_API[0]
ho2=Ho2_stage[23]
quick_plot(ho2)
len(ho2)
ho2=Ho2_stage[23][20:-8]  
'''

stage 6 is 3:-8
stage 22 is 8:-8
stage 23 is 20:-8
'''

Ho2_mean_API[6]=Ho2_stage[6][3:-8].mean()
Ho2_mean_API[22]=Ho2_stage[22][8:-8].mean()
Ho2_mean_API[23]=Ho2_stage[23][20:-8].mean()
Ho2_mean_API_3= Ho2_mean_API
#%% plot the calibration factor for SA
path = os.getcwd()

file=os.path.join(path, "input_files/SA_model_3.csv")
model=pd.read_csv(file)
Sa_model=model['SA'][0:11]
Sa_meas=[Sa_mean_API_3[i] for i in range(3,14)]


f1 = np.polyfit(Sa_meas, Sa_model,1)
f = np.poly1d(f1)
print('%e'%f[0])
print('%e'%f[1])
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Sa_meas, Sa_model)
print('%e'%slope)
print('%e'%r_value**2)


fig,ax=plt.subplots(figsize=(5,4))
plt.style.use('default')
plt.rcParams.update({'font.size':16,'font.weight':'bold','font.family':'serif','font.serif':'Times New Roman'})
ax.scatter(Sa_meas,Sa_model,marker='s',label='SA')
ax.plot(np.unique(Sa_meas),f(np.unique(Sa_meas)),'k',label= r'y=1.6e10*x-2.8e6')

# ax.grid(True)
# ax.grid(which='both')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-4,-5))
ax.yaxis.major.formatter._useMathText = True
# ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_title('SA cali 3')

ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
tkw = dict(size=6, width=1)
ax.tick_params(axis='y', which='major',direction='out',right=True,colors='black', **tkw)
ax.tick_params(axis='x', which='major', direction='out',top=True,colors='black',**tkw)
tkw = dict(size=3, width=1)
ax.tick_params(axis='y', which='minor',direction='out',right=True,colors='black', **tkw)
ax.tick_params(axis='x', which='minor', direction='out',top=True,colors='black',**tkw)

ax.set_xlabel('Norm. $\mathregular{H_2SO_4}$ signal (cps $\mathregular{cps^{-1}}$)',fontweight='bold')
ax.set_ylabel('Modelled $\mathregular{H_2SO_4}$ conc. (molec $\mathregular{cm^{-3}}$)',fontweight='bold')
ax.annotate("Slope = 4.5 × $\mathdefault{10^{9}}$ \n $\mathdefault{R^2}$ = 0.99",
            xy=(.04, .78), xycoords=ax.transAxes)
plt.savefig(path+'/output/SA_model_vs_norma_3.png',bbox_inches='tight', dpi=100)
#%% plot the calibration factor for HO2
path = os.getcwd()

file=os.path.join(path, "input_files/SA_model_3.csv")
model=pd.read_csv(file)
Ho2_model=model['HO2'][:3]
Ho2_meas=[Ho2_mean_API_3[i] for i in range(3,6)]


f1 = np.polyfit(Ho2_meas, Ho2_model,1)
f = np.poly1d(f1)
print('%e'%f[0])
print('%e'%f[1])
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Ho2_meas, Ho2_model)
print('%e'%slope)
print('%e'%r_value**2)


fig,ax=plt.subplots(figsize=(5,4))
plt.style.use('default')
plt.rcParams.update({'font.size':16,'font.weight':'bold','font.family':'serif','font.serif':'Times New Roman'})
ax.scatter(Ho2_meas,Ho2_model,marker='d',label='SA')
ax.plot(np.unique(Ho2_meas),f(np.unique(Ho2_meas)),'k')

# ax.grid(True)
# ax.grid(which='both')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-4,-5))
ax.yaxis.major.formatter._useMathText = True
# ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_title('HO2 cali 3')

ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
tkw = dict(size=6, width=1)
ax.tick_params(axis='y', which='major',direction='out',right=True,colors='black', **tkw)
ax.tick_params(axis='x', which='major', direction='out',top=True,colors='black',**tkw)
tkw = dict(size=3, width=1)
ax.tick_params(axis='y', which='minor',direction='out',right=True,colors='black', **tkw)
ax.tick_params(axis='x', which='minor', direction='out',top=True,colors='black',**tkw)


ax.set_xlabel('Norm. $\mathregular{HO_2}$ signal (cps $\mathregular{cps^{-1}}$)',fontweight='bold')
ax.set_ylabel('Modelled $\mathregular{HO_2}$ conc. (molec $\mathregular{cm^{-3}}$)',fontweight='bold')
ax.annotate("Slope = 3.7 × $\mathdefault{10^{11}}$ \n $\mathdefault{R^2}$ = 0.91, H2O[100,200,400]",
            xy=(.04, .78), xycoords=ax.transAxes)
plt.savefig(path+'/output/HO2_model_vs_norma_3.png',bbox_inches='tight', dpi=100)


#%% plot the RH effect
path = os.getcwd()

file=os.path.join(path, "input_files/SA_model_3.csv")
model=pd.read_csv(file)
Sa_model=model['SA'][10:]
Sa_meas=[Sa_mean_API_3[i]*4.518916e+09 for i in range(13,24)]
h2o_p=[H2O_stage_conc_3[i] for i in range(13,24)] #0,14// 15,22

f1 = np.polyfit(Sa_meas, Sa_model,1)
f = np.poly1d(f1)
print('%e'%f[0])
print('%e'%f[1])
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(Sa_meas, Sa_model)
print('%e'%slope)
print('%e'%r_value**2)


fig,ax=plt.subplots(figsize=(5,4))
plt.style.use('default')
plt.rcParams.update({'font.size':16,'font.weight':'bold','font.family':'serif','font.serif':'Times New Roman'})
ax.scatter(h2o_p,Sa_meas,marker='s',label='SA cali 2')
# ax.plot(np.unique(Sa_meas),f(np.unique(Sa_meas)),'k',label= r'y=1.6e10*x-2.8e6')
# ax.plot()
# ax.grid(True)
# ax.grid(which='both')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-4,-5))
ax.yaxis.major.formatter._useMathText = True
# ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_title('SA RH 3')

ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
tkw = dict(size=6, width=1)
ax.tick_params(axis='y', which='major',direction='out',right=True,colors='black', **tkw)
ax.tick_params(axis='x', which='major', direction='out',top=True,colors='black',**tkw)
tkw = dict(size=3, width=1)
ax.tick_params(axis='y', which='minor',direction='out',right=True,colors='black', **tkw)
ax.tick_params(axis='x', which='minor', direction='out',top=True,colors='black',**tkw)
ax.set_xlim(1e16,4e17)
ax.set_ylim(1e7,4e7)
ax.set_xscale('log')

ax.set_xlabel('$\mathregular{H_2O}$ conc. (molec $\mathregular{cm^{-3}}$)',fontweight='bold')
ax.set_ylabel('$\mathregular{H_2SO_4}$ conc. (molec $\mathregular{cm^{-3}}$)',fontweight='bold')
# ax.annotate("Slope = 7 × $\mathdefault{10^{10}}$ \n $\mathdefault{R^2}$ = 0.95",
#             xy=(.04, .78), xycoords=ax.transAxes)
plt.savefig(path+'/output/SA_vs_H2O_RH_3.png',bbox_inches='tight', dpi=100)

