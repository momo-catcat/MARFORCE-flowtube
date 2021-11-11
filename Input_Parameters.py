%% input parameters
% set these inputs and run the file

addpath('C:\Users\jiali\OneDrive - University of Helsinki\CLOUD15T\calibration\SA_cali_2_RH_Aver\SA_calibrator_Miska\')
% load the water flow from the python data 
H2O=readtable('H2O.csv')

% if there are two different tube then go to the cmd_cali1Matlab.m 
R1= 0.78*4/3; % 1 inch for 58.5 cm 1 inch is 2.54 cm, 0.2 for the tube wall 
L1= 58.5; % cm
R2 = 0.78; % 3/4 inch for 41 cm
L2 = 41; % cm



T0= 273.15;
t= 20; %k
T = T0+t; % K
p = 96060*1.005; % Pa

ID = R1*10; % mm the dimaters of the tube
L = (L2+L1)*10; % mm
Q = 20; % lpm, not slpm

Itx = 5.2009e10; % at Qx flow rate
Qx = 20; % lpm

N2Flow = 22; % slpm
AirFlow = 50; % smlpm
WaterFlow = H2O.H2OCali_MFC4_MFC2_ 
% [0 0 1000 500 750 1000 1250 1500 1750 2000 1500 500 0]'; % smlpm
SO2Flow = 5; % smlpm
SO2BottlePpm = 5000; % ppm

O2inAir = 0.209;

outflowLocation = 'before'; % outflow tube located before or after injecting air, water, and so2

fullOrSimpleModel = 'full'; % simple: Gormley&Kennedy approximation, full: flow model (much slower)

t1=L1/(Q*1e3/60/pi/R1^2);
t2=L2/(Q*1e3/60/pi/R2^2);
%save('SA_model_4.mat','H2SO4','ID','L');%want to save the file for
%different input parameters because we have to differet diameters for the
%inlet
%% computation begins

if strcmp(outflowLocation,'after')
    totFlow = N2Flow+AirFlow/1000+WaterFlow/1000+SO2Flow/1000;
else
    totFlow = Q*ones(size(WaterFlow));
end
O2conc = O2inAir*AirFlow/1000./totFlow*p/1.3806488e-23/T/1e6;
H2Oconc = WaterFlow/1000./totFlow*vappresw(T)/1.3806488e-23/T/1e6;
SO2conc = SO2Flow/1000./totFlow*SO2BottlePpm*1e-6*p/1.3806488e-23/T/1e6;

It = Itx*Qx/Q;

H2SO4 = zeros(size(WaterFlow));
for i=1:numel(H2SO4)
    H2SO4(i)=cmd_calib1Matlab(O2conc(i),H2Oconc(i),SO2conc(i),ID/10,L/10,Q*1000/60,It,T,p,fullOrSimpleModel,t1,t2);
end