%% input parameters
% set these inputs and run the file


T = 301.15; % K
p = 101000; % Pa

ID = 10; % mm
L = 395; % mm
Q = 8.5; % lpm, not slpm

Itx = 6.186e10; % at Qx flow rate
Qx = 7.6; % lpm

N2Flow = 10.5; % slpm
AirFlow = 21.9; % smlpm
WaterFlow = [0.05426566
0.16198066
0.21583816
0.32355316
0.43126816
0.53898316
0.659392
0.768799
0.987613
1.206427] * 1000'; % smlpm
SO2Flow = 10.7; % smlpm
SO2BottlePpm = 1000; % ppm

% 
% %% input parameters
% % set these inputs and run the file
% 
% 
% T = 298.15; % K
% p = 96060*1.005; % Pa
% 
% ID = 7.8*2; % mm
% L = 1500; % mm
% Q = 11; % lpm, not slpm
% 
% Itx = 5.2009e10; % at Qx flow rate
% Qx = 20; % lpm
% 
% N2Flow = 11; % slpm
% AirFlow = 50; % smlpm
% WaterFlow = [2500]'; % smlpm
% SO2Flow = 5; % smlpm
% SO2BottlePpm = 5000; % ppm

O2inAir = 0.209;

outflowLocation = 'before'; % outflow tube located before or after injecting air, water, and so2

fullOrSimpleModel = 'full'; % simple: Gormley&Kennedy approximation, full: flow model (much slower)

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

% H2Oconc = 1.75011e+16;
% O2conc = 1.10846e+16;
% SO2conc = 2.65181e+13;


H2SO4 = zeros(size(WaterFlow));

for i=1:numel(H2SO4)
    [H2SO4(i),c] =cmd_calib1Matlab(O2conc(i),H2Oconc(i),SO2conc(i),ID/10/2,L/10,Q*1000/60,It,T,p,fullOrSimpleModel);
end

disp(H2SO4)

%%
save('c_matlab.mat','c','H2SO4')