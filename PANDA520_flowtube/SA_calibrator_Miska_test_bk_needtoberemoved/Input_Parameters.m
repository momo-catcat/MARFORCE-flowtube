%% input parameters
% set these inputs and run the file


T = 298; % K
p = 101000; % Pa

ID = 7.8*2; % mm
L = 260; % mm
Q = 10.6; % lpm, not slpm

Itx = 5.42e10; % at Qx flow rate
Qx = 20 % lpm

N2Flow = 10.5; % slpm
AirFlow = 50; % smlpm
WaterFlow = [0.100
0.2
0.4
0.8
1
1.2
1.5
0.9
0.6] * 1000'; % smlpm
SO2Flow = 5; % smlpm
SO2BottlePpm = 5000; % ppm

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