import numpy as np

from Vapour_calc import H2O_conc


def const_comp_conc_cal(H2O_data, outflowLocation, const_comp_pre, pre_flow, N2Flow, pre_standard_conc,
                        Q1, Q2, T_cel, T, p, flag_tube):

    kB = 1.3806488e-23  # boltzmann constant


    WaterFlow1 = H2O_data['H2O_1']
    var_flow = [WaterFlow1]

    if outflowLocation in 'after':
        totFlow1 = N2Flow + np.sum(pre_flow) / 1000 + np.sum(var_flow, axis=0) / 1000
    else:
        totFlow1 = Q1

    if 'H2Oconc_1' in H2O_data.columns:
        H2Oconc1 = H2O_data['H2Oconc_1']
    else:
        H2Oconc1 = WaterFlow1 / 1000 / totFlow1 * H2O_conc(T_cel, 1).SatP[0] / kB / T / 1e6

    pre_conc1 = [pre_flow[i] * pre_standard_conc[i] / 1000 / totFlow1 * p / kB / T / 1e6 for i in range(len(pre_flow))]
    var_conc1 = H2Oconc1

    if 'H2O_2' in H2O_data.columns:
        WaterFlow2 = H2O_data['H2O_2']  # Y pieces
    else:
        WaterFlow2 = H2O_data['H2O_1']  # no Y pieces

    totFlow2 = Q2

    if flag_tube in ['1', '2']:
        H2Oconc2 = H2Oconc1
        pre_conc2 = pre_conc1
    else:
        pre_conc2 = [pre_conc1[i] * Q1 / Q2 for i in range(len(pre_conc1))]
        if 'H2Oconc_2' in H2O_data.columns:
            H2Oconc2 = H2O_data['H2Oconc_2']
        else:
            H2Oconc2 = (WaterFlow2 + WaterFlow1) / 1000 / totFlow2 * H2O_conc(T_cel, 1).SatP[0] / kB / T / 1e6

    var_conc2 = H2Oconc2

    var = np.transpose([var_conc1, var_conc2])
    #var = np.asarray(var)
    pre = [pre_conc1, pre_conc2]
    #pre = np.transpose(pre)
    # % store all the const species to const_comp_conc
    #const_comp_conc_values = np.concatenate((pre, var), axis=1)

    if flag_tube == '3':
        const_comp_free = ['H2O', 'O2']
        const_comp_conc_free = [var_conc2[const_comp_var.index('H2O')], pre_conc2[const_comp_pre.index('O2')]]
    else:
        const_comp_free = []
        const_comp_conc_free = [0]

    Itx = 5.2009e10  # at Qx flow rate
    Qx = 20  # lpm
    csH2O = 7.22e-20  # cm2
    qyH2O = 1
    It = Itx * Qx / Q1

    OHconc = It * csH2O * qyH2O * var_conc1[const_comp_var.index('H2O')]

    return var, pre, const_comp_free, const_comp_conc_free, OHconc
