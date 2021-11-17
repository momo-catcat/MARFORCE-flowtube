def odesolve(timesteps, Zgrid, Rgrid, dt, kSO2pOH, kOHpHO2, kOHpOH, kSO3p2H2O, kHSO3pO2, O2conc, H2Oconc, SO2conc, D, R, L, Q, cc):
    #import packages
    import numpy as np
    c = cc
    initc = c

    dz = L / (Zgrid - 1)
    dr = R / (Rgrid - 1)

    for m in range(timesteps):
        term1 = np.zeros([int(Rgrid / 2), int(Zgrid), 5])
        term2 = np.zeros([int(Rgrid / 2), int(Zgrid), 5])
        term3 = np.zeros([int(Rgrid / 2), int(Zgrid), 5])
        for l in range(5):
            for i in range(1, int(Zgrid - 1), 1):
                for j in range(1, int(Rgrid / 2), 1):

                    r = (Rgrid - j) * dr
                    r = np.abs(2 * r - R)

                    term1[j, i, l] = D[l] * (1 / r * (initc[j, i, l] - initc[j - 1, i, l]) / (-2 * dr) + (initc[j + 1, i, l] - 2. * initc[j, i, l] + initc[j - 1, i, l]) / \
                                    (4 * dr * dr) + (initc[j, i + 1, l] - 2. * initc[j, i, l] + initc[j, i - 1, l]) / (dz * dz))

                    term2[j, i, l] = (2 * Q) / (np.pi * R ** 2) * (1. - r ** 2 / (R ** 2)) * (initc[j, i, l] - initc[j, i - 1, l]) / dz

                    if l == 0:
                        term3[j, i, 0] = kSO2pOH * SO2conc * initc[j, i, 4] - kHSO3pO2 * initc[j, i, 0] * O2conc
                    elif l == 1:
                        term3[j, i, 1] = kHSO3pO2 * O2conc * initc[j, i, 0] - kSO3p2H2O * H2Oconc * H2Oconc * initc[j, i, 1]
                    elif l == 2:
                        term3[j, i, 2] = kHSO3pO2 * O2conc * initc[j, i, 0] - kOHpHO2 * initc[j, i, 2] * initc[j, i, 4]
                    elif l == 3:
                        term3[j, i, 3] = kSO3p2H2O * H2Oconc * H2Oconc * initc[j, i,1]
                    elif l == 4:
                        term3[j, i, 4] = -kSO2pOH * SO2conc * initc[j, i, 4] - kOHpHO2 * initc[j, i, 2] * initc[j, i, 4] - 2 * kOHpOH * initc[j, i, 4] * initc[j, i, 4]

            i = Zgrid - 1
            for j in range(1, (Rgrid / 2).astype(int), 1):
                r = (Rgrid - j) * dr
                r = np.abs(2 * r - R)

                term1[j, i, l] = D[l] * (1. / r * (initc[j, i, l] - initc[j - 1, i, l]) / (-2 * dr) + (initc[j + 1, i, l] - 2. * initc[j, i, l] + initc[j - 1, i, l]) / \
                                (4 * dr * dr) + (initc[j, i, l] - 2. * initc[j, i-1, l] + initc[j, i - 2, l]) / (dz * dz))

                term2[j, i, l] = (2. * Q) / (np.pi * R ** 2) * (1. - r ** 2 / (R ** 2)) * (initc[j, i, l] - initc[j, i - 1, l]) / dz

                if l == 0:
                    term3[j, i, 0] = kSO2pOH * SO2conc * initc[j, i, 4] - kHSO3pO2 * initc[j, i, 0] * O2conc
                elif l == 1:
                    term3[j, i, 1] = kHSO3pO2 * O2conc * initc[j, i, 0] - kSO3p2H2O * H2Oconc * H2Oconc * initc[j, i, 1]
                elif l == 2:
                    term3[j, i, 2] = kHSO3pO2 * O2conc * initc[j, i, 0] - kOHpHO2 * initc[j, i, 2] * initc[j, i, 4]
                elif l == 3:
                    term3[j, i, 3] = kSO3p2H2O * H2Oconc * H2Oconc * initc[j, i, 1]
                elif l == 4:
                    term3[j, i, 4] = -kSO2pOH * SO2conc * initc[j, i, 4] - kOHpHO2 * initc[j, i, 2] * initc[
                        j, i, 4] - 2 * kOHpOH * initc[j, i, 4] * initc[j, i, 4]


        c[0 : int(c[:, 0, 0].size / 2), :, :] = dt * (term1 - term2 + term3) + initc[0 : int(c[:, 0, 0].size / 2), :, :]

        c[int(Rgrid / 2) + 1, i, l] = c[int(Rgrid / 2), i, l]


        initc = c

    conc = c
    return(conc)

