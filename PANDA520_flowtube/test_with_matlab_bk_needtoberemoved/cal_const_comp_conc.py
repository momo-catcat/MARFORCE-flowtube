# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 11:00:51 2022

@author: jiali
"""
import numpy as np


## this file is used to set the const comp conc for the grids
def cal_const_comp_conc(Rgrid, Zgrid, const_comp_conc, L1, L2, const_comp):
    const_comp_gird = []
    for i in range(len(const_comp)):
        # const_comp[i]
        H2Otot = np.zeros([int(Rgrid), int(Zgrid)])
        H2Otot[:, 0:int(Zgrid * L1 / (L2 + L1))] = const_comp_conc[0, i]
        H2Otot[:, int(Zgrid * L1 / (L2 + L1)):] = const_comp_conc[1, i]
        const_comp_gird.append(H2Otot)

    return const_comp_gird


## this file is used to set the const comp conc for the grids
def cal_const_comp_conc_1(Rgrid, Zgrid, const_comp_conc, const_comp):
    const_comp_gird = []
    for i in range(len(const_comp)):
        # const_comp[i]
        H2Otot = np.zeros([int(Rgrid), int(Zgrid)])
        H2Otot[:, :] = const_comp_conc[0, i]
        const_comp_gird.append(H2Otot)

    return const_comp_gird
