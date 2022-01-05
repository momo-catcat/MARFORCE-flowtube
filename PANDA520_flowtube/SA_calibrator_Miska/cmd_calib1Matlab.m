function meanWeightedH2SO4=cmd_calib1Matlab(O2conc,H2Oconc,SO2conc,R,L,Q,It,T,p,fullOrSimpleModel) 
    %% parameters (that can be changed by the user)

    if H2Oconc==0
        meanWeightedH2SO4=0;
        return
    end
    % reaction constants (using cm^3 and s)
    kSO2pOH = 1.32e-12*(T/300)^-0.7;
    kOHpHO2 = 4.8e-11*exp(250/T);
    kOHpOH = 6.9e-31*(T/300)^-0.8*p/1.3806488e-23/T/1e6;
    kSO3p2H2O = 3.9e-41*exp(6830.6/T);
    kHSO3pO2 = 1.3e-12*exp(-330/T);
    
    

    % initial concentrations, [1/cm^3]
    % O2conc = 1e17;
    % H2Oconc = 1e16;
    % SO2conc = 1e14;
    % OHconc = 1e9;
    
    % AP 10.7.2013: cmd.m script replaced with cmd_v2.m function
    % in order to perform several model runs with an external script
    % OHconc computed (from H2Oconc)
    
    % Initial OH concentration
    csH2O=7.22e-20; %cm2
    qyH2O=1;
    %It=1.84e10;
    OHconc=It*csH2O*qyH2O*H2Oconc;
    
    % diffusion constants, [cm^2/s]
    % order is: HSO3, SO3, HO2, H2SO4, OH
    % (TODO: check D(H2SO4) RH dependance! we go up to 70% for our calibration)
    rh=H2Oconc*1e6*1.3806488e-23*T/vappresw(T);
    D = [0.126 0.126 0.141 diff_sa_rh(298,rh)*1e4 0.215];

    T0=[300 300 298 298 298];
    D=101325/p*D.*((T.^(3/2))./(T0.^(3/2)));

    if strcmp(fullOrSimpleModel,'simple')
        meanWeightedH2SO4=OHconc*gormleyKennedy(D(4)*1e-4,L/100,Q*60/1000);
        return;
    end

    % other parameters
    %R = .953; % R = .5;                             % radius [cm]
    %L = 15; % L = 25;                             % length [cm]
    %Q = 125; % Q = 167;                            % flow rate [cm^3/s]

    % solving parameters
    dt = 0.00001;                       % timestep [s]
    numLoop = 500;                       % number of times to run, a higher number than 1 will plot intermediate results
    timesteps = 10000;                   % number of timesteps, dt * timesteps * numLoop is time elapsed in the final solution
    Zgrid = 40;                         % number of grid points in tube length direction
    Rgrid = 80;                         % number of grid points in tube radius direction

    % do not change
    while gcd(Rgrid - 1, 10) ~= 1 % make sure we don't get a grid point for r = 0 (would give Inf in calculation)
        Rgrid = Rgrid + 1;
    end

    % initial conditions
    % order in c is: HSO3, SO3, HO2, H2SO4, OH
    c = cell(1, 5);
    c = cellfun(@(x)zeros(Rgrid, Zgrid), c, 'UniformOutput', false);
    c{5}(:, 1) = OHconc;                % set [OH] at z = 0
    c{3}(:, 1) = OHconc;                % set [HO2] at z = 0

    %% check and output parameters to file

    dr = R / (Rgrid - 1) * 2;
%     fId = fopen('par.dat', 'w');
%     fprintf(fId, '%d %d %d %g\n', timesteps, Zgrid, Rgrid, dt);
%     fprintf(fId, '%g %g %g %g %g\n', kSO2pOH, kOHpHO2, kOHpOH, kSO3p2H2O, kHSO3pO2);
%     fprintf(fId, '%g %g %g %g\n', O2conc, H2Oconc, SO2conc, OHconc);
%     fprintf(fId, '%g %g %g %g %g\n', D(1), D(2), D(3), D(4), D(5));
%     fprintf(fId, '%g %g %g\n', R, L, Q);
%     fclose(fId);
%     for i = 1:5
%         a = c{i};
%         eval(sprintf('save c%d.dat a -ascii', i - 1));
%     end

    %% calculate and plot results
    figure; % AP 10.7.2013
    for j = 1:numLoop
%         system('odesolve.exe'); % call C-program

        oldH2SO4 = c{4}(:,end);
        c=odesolveMatlab(timesteps,Zgrid,Rgrid,dt,kSO2pOH,kOHpHO2,kOHpOH,kSO3p2H2O,kHSO3pO2,O2conc,H2Oconc,SO2conc,D,R,L,Q,c);
        t = {'HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH'};
        
        for i = 1:5
%             c{i} = load(sprintf('c%d.dat', i - 1)); % read concentrations calculated by C-program
            c{i}(size(c{i},1)/2+1:end-1,:)=flipud(c{i}(2:size(c{i},1)/2,:)); % symmetry, program only returns values for 0..R
        end
        tim = j * timesteps * dt;
        for i = 1:5 % plot
            subplot(2,3,i);pcolor(linspace(0,L,Zgrid),linspace(-R,R,Rgrid),c{i});shading flat
            %colorbar;
            xlabel('L [cm]');
            ylabel('r [cm]');
            if i == 2
                title({sprintf('t = %f s',tim), t{i}}) % also print time in subplot 2 (top-middle)
            else
                title(t{i});
            end
        end
        drawnow
        newH2SO4 = c{4}(:,end);
        disp(sprintf('t = %f, [H2SO4] change = %g', tim, sum(newH2SO4 - oldH2SO4)));
        
        if j>15&&sum(newH2SO4 - oldH2SO4)/sum(oldH2SO4)<1e-5
            break
        end
        
    end

    %% finally, calculate mean [H2SO4]

    x = R:-dr:0;
    y = c{4}(1:size(c{4},1)/2,end)';
    global splineres % this is to be used by f-function
    splineres = spline(x, y);

    % plot concentration profile
    subplot(2,3,6);
    plot(0:0.01:R,ppval(splineres, 0:0.01:R));
    title('[H2SO4] at end of tube');
    xlabel('r [cm]');
    ylabel('Concentration, [cm^{-3}]');

%     numpoints = 400; % grid to use for integration, larger => greater accuracy
%     xy = zeros(numpoints, numpoints);
%     di = R / numpoints;
%     xmesh = repmat(0:di:R,numpoints + 1,1);
%     ymesh = xmesh';
%     mesh = sqrt(xmesh.^2 + ymesh.^2); % mesh of 0:R x 0:R, value of point is distance to origin
%     conc = ppval(splineres, mesh);
%     conc(mesh > R) = 0; % remove points outside of tube
%     conc1=ones(size(conc));
%     conc1(mesh > R) = 0;
%     concWeighted = conc .* (-sqrt(xmesh.^2 + ymesh.^2).^2 ./ R^2 + 1);
%     conc1Weighted = conc1 .* (-sqrt(xmesh.^2 + ymesh.^2).^2 ./ R^2 + 1);
%     meanH2SO4 = (conc(1,1) + 2 * sum(conc(2:end,1)) + 2 * sum(conc(1,2:end)) + 4 * sum(sum(conc(2:end,2:end)))) * di^2 %/ ((conc1(1,1) + 2 * sum(conc1(2:end,1)) + 2 * sum(conc1(1,2:end)) + 4 * sum(sum(conc1(2:end,2:end)))) * di^2)
%     meanWeightedH2SO4 = (concWeighted(1,1) + 2 * sum(concWeighted(2:end,1)) + 2 * sum(concWeighted(1,2:end)) + 4 * sum(sum(concWeighted(2:end,2:end)))) * di^2% / ((conc1Weighted(1,1) + 2 * sum(conc1Weighted(2:end,1)) + 2 * sum(conc1Weighted(1,2:end)) + 4 * sum(sum(conc1Weighted(2:end,2:end)))) * di^2)


    rVec=0:0.001:R;
    cVec=ppval(splineres, 0:0.001:R);
    meanH2SO4=2*0.001/R^2*sum(cVec.*rVec);
    meanWeightedH2SO4=4*0.001/R^2*sum(cVec.*rVec.*(1-rVec.^2/R^2))

    %{
        maxint = ppval(splineres, 0); % get integration bounds
        splineres = spline(y, x); % disc integration with respect to r-axis
        meanH2SO4 = quadl(@(x)f(x), 0, maxint, 1e3) % integrate
    %}
end
