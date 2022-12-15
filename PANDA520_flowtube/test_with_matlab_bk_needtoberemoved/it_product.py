#%% this script is used to calculate It product for the calibration box
# load the N2 NO and N2O data
import sys
import numpy as np
import os
import pandas as pd

sys.path.append("C:/Users/jiali/PANDA520-flowtube/PANDA520_flowtube/")

file = os.getcwd() + '/input_files/input_data_for_it_product.csv'
data = pd.read_csv(file)
#%%
N2_flow=data['dry N2 af']
N2O_flow=data['N2O']
NO_conc=data['NO']
Temp=20
Total=N2_flow+N2O_flow
#%%
NOx=NO_conc #ppb
## calculation
## parameters
T=Temp+ 273.15 #K
K6=2.0*1e-11*np.exp(130/T) #cm3 mol-1 s-1
K7=7.6*1e-11 #cm3 mol-1 s-1
K8=4.3*1e-11 #cm3 mol-1 s-1
K9=6.0*1e-12 #cm3 mol-1 s-1

phi_N2O=1 #quantum yield
sigma_N2O=1.43*1e-19 #absorption cross section cm2

mixing_N2=N2_flow/(N2O_flow+N2_flow)
mixing_N2O=N2O_flow/(N2O_flow+N2_flow)

#inlet flow CI-API-TOF
flow_inlet=20000 #sccm
#flow_inlet_HOxROx=7450; %sccm

#% calculation of It product
It_N2O=(K6*mixing_N2+(K7+K8+K9)*mixing_N2O)*NOx*1e-9/(2*K7*sigma_N2O*phi_N2O*mixing_N2O ** 2);

#calibrate with the inlet flow
It=It_N2O*(N2_flow+N2O_flow)/flow_inlet
#% It_HOxROx=It_N2O*(N2_flow+N2O_flow)/flow_inlet_HOxROx
#% It=It_N2O



#geometry correction

N2O_conc=96060./1.38e-23/(Temp+273.15)/1e6*mixing_N2O
m=len(mixing_N2O)
R=0.94 #cm
#%%
#for k in range(m):
#    funcart = @(x,y) np.exp(-sigma_N2O*N2O_conc(k)*(np.sqrt(R**2-y**2)+x))
#    ymin = @(x) -np.sqrt(R**2-x**2)
#    ymax = @(x) np.sqrt(R**2-x**2)
#    q[k] = integral2(funcart,-R,R,ymin,ymax)

from scipy import integrate
fx = []
for i in range(m):
    f = lambda y, x: np.exp(-sigma_N2O*N2O_conc[i]*(np.sqrt(R**2-y**2)+x))
    fx1 = integrate.dblquad(f, -R, R, lambda x: -np.sqrt(R**2-x**2), lambda x: np.sqrt(R**2-x**2))
    fx.append(fx1[0])

q=np.transpose(fx)

K = q/R**2/np.pi


#%% It final
It_corr=It/K
# It_corr_HOxROx=It_HOxROx./K';
#mean(It_corr.*Total)
It_final=np.median(It_corr)

# It_final_HOxROx=median(It_corr_HOxROx);

It_final_mean=np.mean(It_final)
It_corr_low=min(It_corr)
It_corr_high=max(It_corr)

#%% calculate the finalized It


#figure(1)
#plt.plot(1:length(It_corr_Bromide),It_corr_Bromide)
#hold on
#plot(1:length(It_corr_HOxROx),It_corr_HOxROx)

#hold off


# %calculate statistics. the data contain beginning stage and ending stage
# %which need to be deleted.
# loca_ex=find(time_final>=datenum(2017,01,01,15,46,00) & time_final<=datenum(2017,01,01,16,47,00));

# boxplot(It_final)

# %calculate the stat
# min_It=min(It_final);
# max_It=max(It_final);
# quantile_It=quantile(It_final,[.25 .5 .75]);
# mean_It=mean(It_final);