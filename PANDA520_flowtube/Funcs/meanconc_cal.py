import numpy as np
from scipy import interpolate
import pandas as pd

def meanconc_cal(R2, Rgrid, plot_spec, comp_namelist, c, model_mode, final_output_method):
    dr_final = R2 / (Rgrid - 1) * 2
    x = np.arange(0, R2, dr_final) + dr_final  # match with y_x

    rVec = np.arange(0, R2, 0.001)

    rVec = np.flip(rVec)

    meanConc = []
    meanweightConc = []
    for i in plot_spec:
        y_x = np.flip(c[:int(Rgrid / 2), -1, comp_namelist.index(i)])  # 'SA'

        splineres1 = interpolate.splrep(x, y_x)

        cVec1 = interpolate.splev(rVec, splineres1)

        y_x = c[int(Rgrid / 2):, -1, comp_namelist.index(i)]  # 'SA'

        splineres2 = interpolate.splrep(x, y_x)

        cVec2 = interpolate.splev(rVec,splineres2)

        conc1 = 0.001 / R2 ** 2 * np.sum(cVec1 * rVec)
        conc2 = 0.001 / R2 ** 2 * np.sum(cVec2 * rVec)

        meanConc.append(conc1 + conc2)

        velocity_profile = 1 - (rVec / R2) ** 2  # Parabolic velocity profile (laminar flow)
        conc1_flow = 2 * 0.001 / (R2 ** 2) * np.sum(cVec1 * rVec * velocity_profile)
        conc2_flow = 2 * 0.001 / (R2 ** 2) * np.sum(cVec2 * rVec * velocity_profile)
        meanweightConc.append(conc1_flow + conc2_flow)

        ####!!!!!!!!!!!!!!!!!Temporary code
        if model_mode == 'kinetic':
            if i == 'H2SO4':
                prof_conc = pd.DataFrame({'R': rVec, 'SA': cVec1})
                prof_conc.to_csv('./Export_files/Theoretical_model.csv')
        ####!!!!!!!!!!!!!!!!!!Temporary code
    if final_output_method == 'mean':
        return meanConc
    elif final_output_method == 'weighted':
        return meanweightConc


def meanconc_cal_sim(R1, Rgrid, comp_namelist, c):
    dr_final = R1 / (Rgrid - 1) * 2
    x = np.arange(0, R1, dr_final) + dr_final
    rVec = np.arange(0, R1, 0.001)

    meanConc = []
    for i in comp_namelist:
        y_x = np.flip(c[0: int(Rgrid / 2), -1, comp_namelist.index(i)])  # 'SA'

        splineres1 = interpolate.splrep(x, y_x)

        cVec = interpolate.splev(rVec, splineres1)
        meanConc.append(2 * 0.001 / R1 ** 2 * np.sum(cVec * rVec))
    return meanConc

