%% input parameters
% set these inputs and run the file


T = 298.15; % K
p = 96060*1.005; % Pa

ID = 24; % mm
L = 930; % mm
Q = 22; % lpm, not slpm

Itx = 5.2009e10; % at Qx flow rate
Qx = 20; % lpm

N2Flow = 22; % slpm
AirFlow = 50; % smlpm
WaterFlow = [50]'; % smlpm
SO2Flow = 5; % smlpm
SO2BottlePpm = 5000; % ppm

O2inAir = 0.209;

outflowLocation = 'after'; % outflow tube located before or after injecting air, water, and so2

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

H2SO4 = zeros(size(WaterFlow));

for i=1:numel(H2SO4)
    H2SO4(i)=cmd_calib1Matlab(O2conc(i),H2Oconc(i),SO2conc(i),ID/10/2,L/10,Q*1000/60,It,T,p,fullOrSimpleModel);
end