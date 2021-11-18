def odesolve(timesteps, Zgrid, Rgrid, dt, kSO2pOH, kOHpHO2, kOHpOH, kSO3p2H2O, kHSO3pO2, O2conc, H2Oconc, SO2conc, D, R, L, Q, cc, OHconc):
    #import packages
    import numpy as np
    import matplotlib.pyplot as plt
    c = cc
    initc = c

    dz = L / (Zgrid - 1)
    dr = R / (Rgrid - 1)
    r = np.zeros([Rgrid, Zgrid, 5])
    r = np.abs(np.array([r[i, :, :] + i - r.shape[0] / 2 + 0.5 for i in range(r.shape[0])])) * 2 * dr

    #do calculations in matrix
    D_temp = np.tile(D, (Rgrid * Zgrid))
    D = np.transpose(np.reshape(D_temp, (Rgrid, Zgrid, 5)), (0, 1, 2))


    #define production and loss terms
    term1 = np.zeros([int(Rgrid), int(Zgrid), 5])
    term2 = np.zeros([int(Rgrid), int(Zgrid), 5])
    term3 = np.zeros([int(Rgrid), int(Zgrid), 5])

    for m in range(timesteps):
        # np.divide(a, b, out=np.zeros_like(a), where=b!=0)

        a = 1. / r[1:-1, 1:-1, :] * (initc[1:-1, 1:-1, :] - initc[0:-2, 1:-1, :]) / (-2 * dr)
        b = (initc[2:, 1:-1, :] - 2. * initc[1:-1, 1:-1, :] + initc[0:-2, 1:-1, :]) / (4 * dr * dr)
        c = (initc[1:-1, 2:, :] - 2. * initc[1:-1, 1:-1, :] + initc[1:-1, 0:-2,: ]) / (dz * dz)

        term1[1:-1, 1:-1, :] = D[1:-1, 1:-1, :] * (a + b + c)

        term2[1:-1, 1:-1, :] = (2. * Q) / (np.pi * R ** 2) * (1. - r[1:-1, 1:-1, :] ** 2 / (R ** 2)) * (initc[1:-1, 1:-1, :] - initc[1:-1, 0:-2, :]) / dz

        # t = ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']
        term3[1:-1, 1:-1, 0] = kSO2pOH * SO2conc * initc[1:-1, 1:-1, 4] - kHSO3pO2 * initc[1:-1, 1:-1, 0] * O2conc

        term3[1:-1, 1:-1, 1] = kHSO3pO2 * O2conc * initc[1:-1, 1:-1, 0] - kSO3p2H2O * H2Oconc * H2Oconc * initc[1:-1, 1:-1, 1]

        term3[1:-1, 1:-1, 2] = kHSO3pO2 * O2conc * initc[1:-1, 1:-1, 0] - kOHpHO2 * initc[1:-1, 1:-1, 2] * initc[1:-1, 1:-1, 4]

        term3[1:-1, 1:-1, 3] = kSO3p2H2O * H2Oconc * H2Oconc * initc[1:-1, 1:-1, 1]

        term3[1:-1, 1:-1, 4] = -kSO2pOH * SO2conc * initc[1:-1, 1:-1, 4] - kOHpHO2 * initc[1:-1, 1:-1, 2] * initc[
            1:-1, 1:-1, 4] - 2. * kOHpOH * initc[1:-1, 1:-1, 4] * initc[1:-1, 1:-1, 4]

        # for l in range(5):
        #     for i in range(1, int(Zgrid - 1), 1):
        #         for j in range(1, int(Rgrid / 2), 1):

                    # r = (Rgrid - j) * dr
                    # r = np.abs(2 * r - R)
                    #
                    # # term1[j, i, l] = D[l] * (1 / r * (initc[j, i, l] - initc[j - 1, i, l]) / (-2 * dr) + (initc[j + 1, i, l] - 2. * initc[j, i, l] + initc[j - 1, i, l]) / \
                    # #                 (4 * dr * dr) + (initc[j, i + 1, l] - 2. * initc[j, i, l] + initc[j, i - 1, l]) / (dz * dz))
                    # #
                    # # term2[j, i, l] = (2 * Q) / (np.pi * R ** 2) * (1. - r ** 2 / (R ** 2)) * (initc[j, i, l] - initc[j, i - 1, l]) / dz
                    #
                    # if l == 0:
                    #     term3[j, i, 0] = kSO2pOH * SO2conc * initc[j, i, 4] - kHSO3pO2 * initc[j, i, 0] * O2conc
                    # elif l == 1:
                    #     term3[j, i, 1] = kHSO3pO2 * O2conc * initc[j, i, 0] - kSO3p2H2O * H2Oconc * H2Oconc * initc[j, i, 1]
                    # elif l == 2:
                    #     term3[j, i, 2] = kHSO3pO2 * O2conc * initc[j, i, 0] - kOHpHO2 * initc[j, i, 2] * initc[j, i, 4]
                    # elif l == 3:
                    #     term3[j, i, 3] = kSO3p2H2O * H2Oconc * H2Oconc * initc[j, i,1]
                    # elif l == 4:
                    #     term3[j, i, 4] = -kSO2pOH * SO2conc * initc[j, i, 4] - kOHpHO2 * initc[j, i, 2] * initc[j, i, 4] - 2 * kOHpOH * initc[j, i, 4] * initc[j, i, 4]
            #
            # i = Zgrid - 1
            # for j in range(1, (Rgrid / 2).astype(int), 1):
            #     r = (Rgrid - j) * dr
            #     r = np.abs(2 * r - R)
            #
            #     term1[j, i, l] = D[l] * (1. / r * (initc[j, i, l] - initc[j - 1, i, l]) / (-2 * dr) + (initc[j + 1, i, l] - 2. * initc[j, i, l] + initc[j - 1, i, l]) / \
            #                     (4 * dr * dr) + (initc[j, i, l] - 2. * initc[j, i-1, l] + initc[j, i - 2, l]) / (dz * dz))
            #
            #     term2[j, i, l] = (2. * Q) / (np.pi * R ** 2) * (1. - r ** 2 / (R ** 2)) * (initc[j, i, l] - initc[j, i - 1, l]) / dz
            #
            #     if l == 0:
            #         term3[j, i, 0] = kSO2pOH * SO2conc * initc[j, i, 4] - kHSO3pO2 * initc[j, i, 0] * O2conc
            #     elif l == 1:
            #         term3[j, i, 1] = kHSO3pO2 * O2conc * initc[j, i, 0] - kSO3p2H2O * H2Oconc * H2Oconc * initc[j, i, 1]
            #     elif l == 2:
            #         term3[j, i, 2] = kHSO3pO2 * O2conc * initc[j, i, 0] - kOHpHO2 * initc[j, i, 2] * initc[j, i, 4]
            #     elif l == 3:
            #         term3[j, i, 3] = kSO3p2H2O * H2Oconc * H2Oconc * initc[j, i, 1]
            #     elif l == 4:
            #         term3[j, i, 4] = -kSO2pOH * SO2conc * initc[j, i, 4] - kOHpHO2 * initc[j, i, 2] * initc[
            #             j, i, 4] - 2 * kOHpOH * initc[j, i, 4] * initc[j, i, 4]




        # c = dt * (term1 - term2 + term3) + initc

        c = dt * (term1 - term2 + term3) + initc
        c[np.abs(c < 1e-30)] = 0
        c[np.abs(c > 1e100)] = 0

        initc = c

    conc = c
    return(conc)

