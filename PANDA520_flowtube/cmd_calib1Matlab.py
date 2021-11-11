def meanWeightedH2SO4=cmd_calib1Matlab(O2conc,H2Oconc,SO2conc,R,L,Q,It,T,p,fullOrSimpleModel,t1,t2):
    ## parameters (that can be changed by the user)

    if H2Oconc == 0:
        meanWeightedH2SO4=0

    # reaction constants (using cm**3 and s)
    kSO2pOH = 1.32e-12 * (T / 300) ** -0.7
    kOHpHO2 = 4.8e-11 * np.exp(250 / T)
    kOHpOH = 6.9e-31 * (T / 300) ** -0.8 * p / 1.3806488e-23 / T / 1e6
    kSO3p2H2O = 3.9e-41 * np.exp(6830.6 / T)
    kHSO3pO2 = 1.3e-12* np.exp(-330 / T)
    
    # Initial OH concentration
    csH2O = 7.22e-20 #cm2
    qyH2O = 1
    #It=1.84e10
    OHconc = It * csH2O * qyH2O * H2Oconc
    
    # diffusion constants, [cm**2/s]
    # order is: HSO3, SO3, HO2, H2SO4, OH
    rh = H2Oconc * 1e6 * 1.3806488e-23 * T / vappresw(T)
    D = [0.126, 0.126, 0.141, diff_sa_rh(298,rh) * 1e4, 0.215]
    
    T0 = [300 300 298 298 298]
    D = 101325 / p * D * ((T ** (3/2)) / (T0 ** (3/2)))


    dt = 0.00001                       # timestep [s]
    timesteps = 4000                   # number of timesteps, dt * timesteps * numLoop is time elapsed in the final solution
    numLoop = np.floor((t1+t2)/dt/timesteps)                       # number of times to run to reach the pinhole of the instrument

    Zgrid = 40                         # number of grid points in tube length direction
    Rgrid = 80                         # number of grid points in tube radius direction

    # do not change
    while np.gcd(Rgrid - 1, 10) != 1  # make sure we don't get a grid point for r = 0 (would give Inf in calculation)
        Rgrid = Rgrid + 1

    # initial conditions
    # order in c is: HSO3, SO3, HO2, H2SO4, OH
    c = cell(1, 5)
    c = cellfun(@(x)zeros(Rgrid, Zgrid), c, 'UniformOutput', false)
    c{5}(:, 1) = OHconc                # set [OH] at z = 0
    c{3}(:, 1) = OHconc                # set [HO2] at z = 0

    ## check and output parameters to file

    dr = R / (Rgrid - 1) * 2

    plt.figure(1, figsize = [16,9])
    for j = 1:numLoop1
#         system('odesolve.exe') # call C-program

        oldH2SO4 = c{4}(:,end)

        c=odesolveMatlab(timesteps,Zgrid,Rgrid,dt,kSO2pOH,kOHpHO2,kOHpOH,kSO3p2H2O,kHSO3pO2,O2conc,H2Oconc,SO2conc,D,R,L,Q,c)

        t = {'HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH'}

        for i = 1:5
#             c{i} = load(sprintf('c#d.dat', i - 1)) # read concentrations calculated by C-program
            c{i}(size(c{i},1)/2+1:end-1,:)=flipud(c{i}(2:size(c{i},1)/2,:)) # symmetry, program only returns values for 0..R
        end
        tim = j * timesteps * dt
        for i = 1:5 # plot
            subplot(2,3,i)pcolor(linspace(0,L,Zgrid),linspace(-R,R,Rgrid),c{i})shading flat
            #colorbar
            xlabel('L [cm]')
            ylabel('r [cm]')
            if i == 2
                title({sprintf('t = #f s',tim), t{i}}) # also print time in subplot 2 (top-middle)
            else
                title(t{i})
            end
        end
        drawnow
        newH2SO4 = c{4}(:,end)
        disp(sprintf('t = #f, [H2SO4] change = #g', tim, sum(newH2SO4 - oldH2SO4)))


    dr = R / (Rgrid - 1) * 2
    x = R:-dr:0
    y = c{4}(1:size(c{4},1)/2,end)'
    global splineres # this is to be used by f-function
    splineres = spline(x, y)

    # plot concentration profile
    subplot(2,3,6)
    plot(0:0.01:R,ppval(splineres, 0:0.01:R))
    title('[H2SO4] at end of tube')
    xlabel('r [cm]')
    ylabel('Concentration, [cm**{-3}]')


    # rVec=0:0.001:R
    # cVec=ppval(splineres, 0:0.001:R)
    # meanH2SO4=2*0.001/R**2*sum(cVec.*rVec)
    # meanWeightedH2SO4=4*0.001/R**2*sum(cVec.*rVec.*(1-rVec.**2/R**2))
    #
    # #{
    #     maxint = ppval(splineres, 0) # get integration bounds
    #     splineres = spline(y, x) # disc integration with respect to r-axis
    #     meanH2SO4 = quadl(@(x)f(x), 0, maxint, 1e3) # integrate

