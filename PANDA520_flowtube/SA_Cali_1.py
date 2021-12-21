# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 17:36:45 2021

@author: jiali
"""

#%% analysis the data of api_data_1
#directly load the csv of SA_cali_1 file
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

file=os.path.join(path, "input_files/api_data_1.csv")

api_data=pd.read_csv(file)
api_data['time']=pd.to_datetime(api_data['time'])
# api_data['time']=api_data['time']- datetime.timedelta(hours=1)
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

file_name= os.path.join(path, "input_files/SA_cali_1.xlsx")

time_table=pd.read_excel(file_name, index_col=None, header=2,engine='openpyxl')

time_table['Start time']=pd.to_datetime(time_table['Start time'])
time_table['End time']=pd.to_datetime(time_table['End time'])
# check the table remove nan rows
print(time_table)

# detele rows that are not usefull
# time_table=time_table.drop(index=[24,25,26,27]) 

print(time_table.columns)



# check the H2O and UV lights 
H2o=time_table[['H20/N2 (mlpm)','UVC.1 filter']]
# save the water flow to the folder
H2o.to_csv(path+"/input_files/H2O_1.csv")

Sa_stage=[]
H2O_stage=[]
Ho2_stage=[]
for i in range(len(time_table)):
    Sa_stage.append(Sa[time_table['Start time'][i]:time_table['End time'][i]])
    Ho2_stage.append(HO2[time_table['Start time'][i]:time_table['End time'][i]])
    H2O_stage.append(H2O[time_table['Start time'][i]:time_table['End time'][i]])

Sa_mean_API=[Sa_stage[i][10:-10].mean() for i in range(len(Sa_stage))] 
H2O_mean_API=[H2O_stage[i][10:-10].mean() for i in range(len(H2O_stage))] 
Ho2_mean_API=[Ho2_stage[i][10:-10].mean() for i in range(len(Ho2_stage))] 
# det_limit_sa=get_detect_llimit(Sa_stage[0])
#% check the SA measured mean values
quick_plot(Sa)
sa=pd.concat([Sa_stage[1],Sa_stage[2]],axis=0)
Sa_mean_API[11]
sa=Sa_stage[4]
quick_plot(sa)
len(sa)
sa=Sa_stage[4][10:-10]  
quick_plot(sa)
'''
stage 3 is 50:-10
stage 7 is 11:-10
stage 12 is 70:
stage 13 is 30:-10  
stage 15 is 30:-10


'''
Sa_mean_API[3]=Sa_stage[3][50:-10].mean()
Sa_mean_API[4]=Sa_stage[4][11:-10].mean()
Sa_mean_API[12]=Sa_stage[12][70:].mean()
Sa_mean_API[13]=Sa_stage[13][30:-10].mean()
Sa_mean_API[15]=Sa_stage[15][30:-10].mean()

Sa_mean_API_1=Sa_mean_API
#% check the HO2 measured mean values
quick_plot(HO2)
sa=pd.concat([Ho2_stage[1],Ho2_stage[2]],axis=0)
Ho2_mean_API[1]
ho2=Ho2_stage[15]
quick_plot( ho2)
len(HO2)
ho2=Ho2_stage[15][20:-10]  
quick_plot(ho2)
'''
stage 3 is 50:-10
stage 15 is 20:-10


'''
Ho2_mean_API[3]=Ho2_stage[3][50:-10].mean()
Ho2_mean_API[15]=Ho2_stage[15][20:-10].mean()


Ho2_mean_API_1=Ho2_mean_API
#% check the H2O
quick_plot(H2O)
h2o=pd.concat([H2O_stage[1],H2O_stage[2]],axis=0)
H2O_mean_API[13]
h2o=H2O_stage[15]
quick_plot(h2o)
len(h2o)
h2o=H2O_stage[15][15:-10]
h2o.mean()  
quick_plot(h2o)

'''
stage 0 is 25:-1
stage 2 is :-25
stage 3 is 50:-10
stage 7 is 11:-10
stage 15 is 15:-10 


'''
H2O_mean_API[0]=H2O_stage[0][25:-1].mean()
H2O_mean_API[2]=H2O_stage[2][:-25].mean()
H2O_mean_API[3]=H2O_stage[3][50:-10].mean()
H2O_mean_API[7]=H2O_stage[7][11:-10].mean()
H2O_mean_API[15]=H2O_stage[15][15:-10].mean()
H2O_stage_conc_1 =[H2O_mean_API[i] *10 ** 18.02663192 for i in range(len(H2O_mean_API))]
#%% plot the calibration factor for SA
path = os.getcwd()

file=os.path.join(path, "input_files/SA_model_1.csv")
model=pd.read_csv(file)
Sa_model=model['SA']
Sa_meas=[Sa_mean_API_1[i] for i in range(4,15)]


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

ax.grid(True)
ax.grid(which='both')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-4,-5))
ax.yaxis.major.formatter._useMathText = True
# ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_title('SA cali 1')
ax.set_xlabel('Norm. $\mathregular{H_2SO_4}$ signal (cps $\mathregular{cps^{-1}}$)',fontweight='bold')
ax.set_ylabel('Modelled $\mathregular{H_2SO_4}$ conc. (molec $\mathregular{cm^{-3}}$)',fontweight='bold')
ax.annotate("Slope = 3.0 × $\mathdefault{10^{10}}$ \n $\mathdefault{R^2}$ = 0.95",
            xy=(.04, .78), xycoords=ax.transAxes)
plt.savefig(path+'/output/SA_model_vs_norma_1.png',bbox_inches='tight', dpi=100)

#%% plot the calibration factor for HO2
path = os.getcwd()

file=os.path.join(path, "input_files/SA_model_1.csv")
model=pd.read_csv(file)
Ho2_model=model['HO2'][[0,1,2,3,4,5,6,9,10]]
Ho2_meas=[Ho2_mean_API_1[i] for i in [4,5,6,7,8,9,10,13,14]]


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

ax.grid(True)
ax.grid(which='both')
ax.ticklabel_format(axis='y', style='sci', scilimits=(-4,-5))
ax.yaxis.major.formatter._useMathText = True
# ax.set_yscale('log')
# ax.set_xscale('log')
ax.set_title('HO2 cali 1')
ax.set_xlabel('Norm. $\mathregular{HO_2}$ signal (cps $\mathregular{cps^{-1}}$)',fontweight='bold')
ax.set_ylabel('Modelled $\mathregular{HO_2}$ conc. (molec $\mathregular{cm^{-3}}$)',fontweight='bold')
ax.annotate("Slope = 2.0 × $\mathdefault{10^{12}}$ \n $\mathdefault{R^2}$ = 0.92, H2O=[1,10,5,7.5,\n 10,12.5,15,15,5]",
            xy=(.04, .7), xycoords=ax.transAxes)
plt.savefig(path+'/output/HO2_model_vs_norma_1.png',bbox_inches='tight', dpi=100)

#%% plot for SA vs H2O

h2o_p=[H2O_stage_conc_1[i] for i in range(2,16)] #0,13// 13,24
sa_p=[Sa_mean_API_1[i] for i in range(2,16)]
fig,ax=plt.subplots(figsize=(5,4))
plt.style.use('default')
plt.rcParams.update({'font.size':16,'font.weight':'bold','font.family':'serif','font.serif':'Times New Roman'})
ax.scatter(h2o_p,sa_p,label='SA')
ax.grid(True)
ax.grid(which='both')
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylabel('SA signal (API)')
ax.set_xlabel('$\mathregular{H_2O}$ conc')
# ax.xaxis.set_major_locator(ticker.AutoLocator())
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %H:%M'))
plt.savefig(path+'/output/SA_CALI_1_SA_VS_H2O.png')