import numpy as np

from Vapour_calc import H2O_conc


def const_comp_conc_cal(O2flow, SO2flow, outflowLocation, sampflow, H2O_1, N2Flow, O2_ratio, SO2_ratio,
                        Q1, Q2, T_cel, T, p, flag_tube):

    kB = 1.3806488e-23  # boltzmann constant
    sumflow = O2flow + SO2flow + H2O_1

    if outflowLocation in 'after':
        totFlow1 = N2Flow + sumflow / 1000
    elif flag_tube in ['3']:
        totFlow1 = Q1
    else:
        totFlow1 = sampflow

    O2conc1 = O2flow * O2_ratio / 1000 / totFlow1 * p / kB / T / 1e6
    SO2conc1 = SO2flow * SO2_ratio / 1000 / totFlow1 * p / kB / T / 1e6

    if flag_tube in ['1', '2']:
        O2conc2 = O2conc1
        SO2conc2 = SO2conc1
    else:
        O2conc2 = O2conc1 * Q1 / Q2 ####!!!!! The situation when Q2 has oxygen is not considered here
        SO2conc2 = SO2conc1 * Q1 / Q2 ####!!!!! The situation when Q2 has SO2 is not considered here

    return np.transpose([O2conc1, O2conc2]), np.transpose([SO2conc1, SO2conc2])


def const_comp_conc_cal_H2O(O2flow, SO2flow, outflowLocation, sampflow, H2O_1, H2O_2, N2Flow, O2_ratio,
                        Q1, Q2, T_cel, T, p, flag_tube):
    kB = 1.3806488e-23  # boltzmann constant
    sumflow = O2flow + SO2flow + H2O_1

    if outflowLocation in 'after':
        totFlow1 = N2Flow + sumflow / 1000
    elif flag_tube in ['3']:
        totFlow1 = Q1
    else:
        totFlow1 = sampflow

    H2Oconc1 = H2O_1 / 1000 / totFlow1 * H2O_conc(T_cel, 1).SatP[0] / kB / T / 1e6
    totFlow2 = Q2
    if flag_tube in ['1', '2']:
        H2Oconc2 = H2Oconc1
    else:
        H2Oconc2 = (H2O_1 + H2O_2) / 1000 / totFlow2 * H2O_conc(T_cel, 1).SatP[0] / kB / T / 1e6

    return np.transpose([H2Oconc1, H2Oconc2])


def const_comp_conc_cal_OH(H2Oconc, O2conc, Q1, flag_tube, Itx, Qx):

    csH2O = 7.22e-20  # cm2
    qyH2O = 1
    It = Itx * Qx / Q1
    OHconc = It * csH2O * qyH2O * H2Oconc[:, 0]
    if flag_tube == '3':
        const_comp_free = ['H2O', 'O2']
        const_comp_conc_free = [H2Oconc[:, 0], O2conc[:, 0]]
    else:
        const_comp_free = []
        const_comp_conc_free = [0]

    return OHconc, const_comp_free, const_comp_conc_free