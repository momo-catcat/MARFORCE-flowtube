# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 11:00:51 2022

@author: jiali
"""
import numpy as np


## this file is used to set the const comp conc for the grids
def cal_const_comp_conc(const_comp_conc, modelparams):
    modelparams.const_comp_gird = []
    for i in range(len(modelparams.const_comp)):
        # const_comp[i]
        H2Otot = np.zeros([int(modelparams.Rgrid), int(modelparams.Zgrid)])
        H2Otot[:, 0:int(modelparams.Zgrid * modelparams.L1 / (modelparams.L2 + modelparams.L1))] = const_comp_conc[0, i]
        H2Otot[:, int(modelparams.Zgrid * modelparams.L1 / (modelparams.L2 + modelparams.L1)):] = const_comp_conc[1, i]
        modelparams.const_comp_gird.append(H2Otot)

    return modelparams


def set_const_comp_conc(c,modelparams):
    for i in modelparams.const_comp:  # set the constant concentrations for const_comp
        try:
            c[:,:,modelparams.comp_namelist.index(i)] = modelparams.const_comp_gird[modelparams.const_comp.index(i)]
        except:
            pass
    return c


def set_const_comp_conc_for_1st_tube(c,modelparams):
    for i in modelparams.const_comp:  # set the constant concentrations for const_comp
        try:
            c[:,:,modelparams.comp_namelist.index(i)] =modelparams.stage_const_comp_conc[0, modelparams.const_comp.index(i)]
        except:
            pass
    return c

def set_const_comp_conc_for_2nd_tube(c,modelparams):
    for i in modelparams.const_comp:  # set the constant concentrations for const_comp
        try:
            c[:,:,modelparams.comp_namelist.index(i)] =modelparams.stage_const_comp_conc[1, modelparams.const_comp.index(i)]
        except:
            c[:,:,modelparams.comp_namelist.index(i)] =modelparams.stage_const_comp_conc[0, modelparams.const_comp.index(i)]
    return c

