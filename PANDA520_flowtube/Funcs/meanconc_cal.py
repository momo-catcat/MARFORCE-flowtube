import numpy as np
from scipy import interpolate
import pandas as pd

def meanconc_cal(c, modelparams, R=None):
    """Outlet concentration of every species: 'mean' (area average, default) or 'weighted'
    (flow-weighted) according to modelparams.final_output_method."""
    if R is None:
        R = modelparams.R2
    dr_final = R / (modelparams.Rgrid - 1) * 2
    x = np.arange(0, R, dr_final) + dr_final  # match with y_x

    rVec = np.arange(0, R, 0.001)

    rVec = np.flip(rVec)

    meanConc = []
    for i in modelparams.comp_namelist:
        y_x = np.flip(c[:int(modelparams.Rgrid / 2), -1, modelparams.comp_namelist.index(i)])  # 'SA'

        splineres1 = interpolate.splrep(x, y_x)

        cVec1 = interpolate.splev(rVec, splineres1)

        y_x = c[int(modelparams.Rgrid / 2):, -1, modelparams.comp_namelist.index(i)]  # 'SA'

        splineres2 = interpolate.splrep(x, y_x)

        cVec2 = interpolate.splev(rVec, splineres2)

        if getattr(modelparams, 'final_output_method', 'mean') == 'weighted':
            # flow-weighted (flux) mean with the parabolic laminar velocity profile, as in the Matlab model
            velocity_profile = 1 - (rVec / R) ** 2
            conc1 = 2 * 0.001 / R ** 2 * np.sum(cVec1 * rVec * velocity_profile)
            conc2 = 2 * 0.001 / R ** 2 * np.sum(cVec2 * rVec * velocity_profile)
        else:
            # 'mean': area-weighted average over the outlet cross-section
            conc1 = 0.001 / R ** 2 * np.sum(cVec1 * rVec)
            conc2 = 0.001 / R ** 2 * np.sum(cVec2 * rVec)

        meanConc.append(conc1 + conc2)

        ####!!!!!!!!!!!!!!!!!Temporary code
        if modelparams.model_mode == 'kinetic':
            if i == 'H2SO4':
                prof_conc = pd.DataFrame({'R': rVec, 'SA': cVec1})
                prof_conc.to_csv('./Export_files/Theoretical_model.csv')
        ####!!!!!!!!!!!!!!!!!!Temporary code

    return meanConc

def meanconc_cal_one_sepcies(c, name, modelparams, R=None):
    if R is None:
        R = modelparams.R2
    dr_final = R / (modelparams.Rgrid - 1) * 2
    x = np.arange(0, R, dr_final) + dr_final  # match with y_x

    rVec = np.arange(0, R, 0.001)

    rVec = np.flip(rVec)

    for i in name:
        y_x = np.flip(c[:int(modelparams.Rgrid / 2), -1, modelparams.comp_namelist.index(i)])  # 'SA'

        splineres1 = interpolate.splrep(x, y_x)

        cVec1 = interpolate.splev(rVec, splineres1)

        y_x = c[int(modelparams.Rgrid / 2):, -1, modelparams.comp_namelist.index(i)]  # 'SA'

        splineres2 = interpolate.splrep(x, y_x)

        cVec2 = interpolate.splev(rVec, splineres2)

        conc1 = 0.001 / R ** 2 * np.sum(cVec1 * rVec)
        conc2 = 0.001 / R ** 2 * np.sum(cVec2 * rVec)

        meanConc = conc1 + conc2
    return meanConc


def meanconc_cal_sim(c,R1,modelparams):
    dr_final = R1 / (modelparams.Rgrid - 1) * 2
    x = np.arange(0, R1, dr_final) + dr_final
    rVec = np.arange(0, R1, 0.001)

    meanConc = []
    for i in modelparams.comp_namelist:
        y_x = np.flip(c[0: int(modelparams.Rgrid / 2), -1, modelparams.comp_namelist.index(i)])  # 'SA'
        splineres1 = interpolate.splrep(x, y_x)

        cVec = interpolate.splev(rVec, splineres1)
        meanConc.append(2 * 0.001 / R1 ** 2 * np.sum(cVec * rVec))
    return meanConc
