def cmd_calib1(O2conc, H2Oconc, SO2conc, R, L, Q, It, T, p, fullOrSimpleModel, time):
    ## parameters (that can be changed by the user)
    #import packages
    import numpy as np
    from Vapour_calc import H2O_conc
    import matplotlib.pyplot as plt
    from odesolve import odesolve as odesolve
    from scipy import interpolate

    import matplotlib
    # matplotlib.use("TkAgg")

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
    D = np.array([0.126, 0.126, 0.141, 0.08, 0.215]) #todo need to calculate Diffusion coefficient properly
    # D = np.array([0.4, 0.4, 0.4, 0.4, 0.4])
    T0 = np.array([300, 300, 298, 298, 298])
    D = 101325 / p * D * ((T ** (3 / 2)) / (T0 ** (3 / 2)))


    dt = 0.00001                        # timestep [s]
    numLoop = 500                      # number of times to run to reach the pinhole of the instrument
    timesteps = 10000                   # number of timesteps, dt * timesteps * numLoop is time elapsed in the final solution

    Zgrid = np.array(40).astype(int)                         # number of grid points in tube length direction
    Rgrid = np.array(80)                         # number of grid points in tube radius direction

    #Change odd number Rgrid to even number grid
    if (Rgrid % 2) != 0:
        Rgrid = Rgrid + 1

    # # do not change
    # while np.gcd(Rgrid - 1, 10) != 1:  # make sure we don't get a grid point for r = 0 (would give Inf in calculation)
    #     Rgrid = int(Rgrid + 1)


    # initial conditions
    # order in c is: HSO3, SO3, HO2, H2SO4, OH
    c = np.zeros([Rgrid, Zgrid, 5])
    c[:, 0, 4] = OHconc  # set [OH] at z = 0
    c[:, 0, 2] = OHconc  # set [HO2] at z = 0. This equals OH conc
    ## check and output parameters to file


    # plt.figure(1, figsize = [8,6])
    for j in range(numLoop):
        oldH2SO4 = c[:, -1 , 3]
        c = odesolve(timesteps, Zgrid, Rgrid, dt, kSO2pOH, kOHpHO2, kOHpOH, kSO3p2H2O, kHSO3pO2, O2conc, H2Oconc, SO2conc, D, R, L, Q, c)

        print(['OH conc1: ' + str(c[40,0,4]) + 'OH conc2: ' + str(c[40,1,4]) + 'OH conc5: ' + str(c[40,10,4])])
        t = ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']

        # for i in range(5):
        #     c[int(c[:, 0, i].size / 2) : -1, :, i] = np.flipud(c[1 : int(c[:, 0, i].size / 2), :, i]) # symmetry, program only returns values for 0..R

        tim = j * timesteps * dt
        print(['time: ' + str(tim)])

        for i in range(5): # plot
            plt.subplot(2, 3, i + 1)
            plt.pcolor(np.linspace(0, L, Zgrid), np.linspace(-R, R, Rgrid), c[:, : ,i], shading = 'nearest')
            # plt.pcolor(c[:, : , i])
            #colorbar
            plt.xlabel('L [cm]')
            plt.ylabel('r [cm]')
            if i == 2:
                plt.title(['t =' + str(tim) + t[i]]) # also print time in subplot 2 (top-middle)
            else:
                plt.title(t[i])
        plt.draw()
        plt.pause(1)


        newH2SO4 = c[:, -1, 3]
        print(['t = ' + str(tim) + "  H2SO4 difference: " + str(np.sum(newH2SO4 - oldH2SO4))])

        if (j > 15) & (np.sum(newH2SO4 - oldH2SO4) / np.sum(oldH2SO4) < 1e-5):
            break

        dr = R / (Rgrid - 1) * 2
        x = np.arange(0, R, dr)
        y = c[0 : int(Rgrid / 2), -1, 3]
        splineres = interpolate.splrep(x, y)

        # plot concentration profile
        plt.subplot(2,3,6)
        plt.plot(np.arange(0, R, dr), interpolate.splev(np.arange(0, R, dr), splineres))
        plt.title('[H2SO4] at end of tube')
        plt.xlabel('r [cm]')
        plt.ylabel('Concentration, [cm**{-3}]')


        rVec = np.arange(0, R, 0.001)
        cVec = interpolate.splev(rVec, splineres)
        meanH2SO4 = 2 * 0.001 / R ** 2 * np.sum(cVec * rVec) #calculate average concentration over the cross section
        meanWeightedH2SO4 = 4 * 0.001 / R ** 2 * np.sum(cVec * rVec * (1 - rVec ** 2 / R ** 2)) #the inner layer shoul have large chance to get into the pinhole
        print(np.sum(cVec * rVec))
        print(np. sum(cVec * rVec * (1 - rVec ** 2 / R ** 2)))

    return(meanWeightedH2SO4)
