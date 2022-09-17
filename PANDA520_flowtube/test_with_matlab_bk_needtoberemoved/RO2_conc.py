# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:45:20 2021

@author: jiali
"""

'''
get the RO2 concentration 
'''


def RO2_conc(RO2_indices, y):
    RO2conc = 0
    for i in range(len(RO2_indices)):
        RO2conc += y[RO2_indices[i, 1]]
    return RO2conc
