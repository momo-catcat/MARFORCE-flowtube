def odesolve(timesteps, Zgrid, Rgrid, dt, kSO2pOH, kOHpHO2, kOHpOH, kSO3p2H2O, kHSO3pO2, O2conc, H2Oconc, SO2conc, D, R, L, Q, c):
    #import packages
    import numpy as np
    import matplotlib.pyplot as plt
    initc = c

    dx = L / (Zgrid - 1)
    dr = 2 * R / (Rgrid - 1)
    r = np.zeros([Rgrid, Zgrid, 5])
    r = np.abs(np.array([r[i, :, :] + i - r.shape[0] / 2 + 0.5 for i in range(r.shape[0])])) * dr

    #do calculations in matrix
    D_temp = np.tile(D, (Rgrid * Zgrid))
    D = np.transpose(np.reshape(D_temp, (Rgrid, Zgrid, 5)), (0, 1, 2))


    #define production and loss terms
    term1 = np.zeros([int(Rgrid), int(Zgrid), 5])
    term2 = np.zeros([int(Rgrid), int(Zgrid), 5])
    term3 = np.zeros([int(Rgrid), int(Zgrid), 5])

    for m in range(timesteps):

        # The equation is based on Fick's first law: https://en.wikipedia.org/wiki/Fick%27s_laws_of_diffusion
        # calculate central grids
        # p_a = 1. / r[1:-1, 1:-1, :] * (initc[1:-1, 1:-1, :] - initc[0:-2, 1:-1, :]) / (-dr). # This doesn't seem to be necessary since we are considering two dimensional problem.
        p_b = (initc[2:, 1:-1, :] - 2. * initc[1:-1, 1:-1, :] + initc[0:-2, 1:-1, :]) / (dr * dr)
        p_c = (initc[1:-1, 2:, :] - 2. * initc[1:-1, 1:-1, :] + initc[1:-1, 0:-2,: ]) / (dx * dx)

        term1[1:-1, 1:-1, :] = D[1:-1, 1:-1, :] * (p_b + p_c) #diffusion of gas molecules

        term2[1:-1, 1:-1, :] = (2. * Q) / (np.pi * R ** 2) * (1. - r[1:-1, 1:-1, :] ** 2 / (R ** 2)) *\
                               (initc[1:-1, 1:-1, :] - initc[1:-1, 0:-2, :]) / dx #carried by main flow

        #calculate the last column (measured by the instrument)
        p_b_end = (initc[2:, -1, :] - 2. * initc[1:-1, -1, :] + initc[0:-2, -1, :]) / (dr * dr)
        p_c_end = (initc[1:-1, -1, :] - 2. * initc[1:-1, -2, :] + initc[1:-1, -2,: ]) / (dx * dx)

        term1[1:-1, -1, :] = D[1:-1, -1, :] * (p_b_end + p_c_end) #diffusion of gas molecules

        term2[1:-1, -1, :] = (2. * Q) / (np.pi * R ** 2) * (1. - r[1:-1, -1, :] ** 2 / (R ** 2)) *\
                               (initc[1:-1, -1, :] - initc[1:-1, -2, :]) / dx #carried by main flow

        # t = ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']
        term3[1:-1, 1:, 0] = kSO2pOH * SO2conc * initc[1:-1, 1:, 4] - kHSO3pO2 * initc[1:-1, 1:, 0] * O2conc

        term3[1:-1, 1:, 1] = kHSO3pO2 * O2conc * initc[1:-1, 1:, 0] - kSO3p2H2O * H2Oconc * H2Oconc * initc[1:-1, 1:, 1]

        term3[1:-1, 1:, 2] = kHSO3pO2 * O2conc * initc[1:-1, 1:, 0] - kOHpHO2 * initc[1:-1, 1:, 2] * initc[1:-1, 1:, 4]

        term3[1:-1, 1:, 3] = kSO3p2H2O * H2Oconc * H2Oconc * initc[1:-1, 1:, 1]

        term3[1:-1, 1:, 4] = -kSO2pOH * SO2conc * initc[1:-1, 1:, 4] - kOHpHO2 * initc[1:-1, 1:, 2] * initc[
            1:-1, 1:, 4] - 2. * kOHpOH * initc[1:-1, 1:, 4] * initc[1:-1, 1:, 4]

        c = dt * (term1 - term2 + term3) + initc

        initc = c

    return(c)

