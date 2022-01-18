# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 09:45:38 2022

@author: jiali
"""

#% help judge the species for reaction rates 
 
def jude_species(y,comp_namelist):
    if 'O2' in comp_namelist:
        O2conc = y[comp_namelist.index('O2')]
    else:
        O2conc = 0
    
    if 'NO' in comp_namelist:
        NOconc = y[comp_namelist.index('NO')]
    else:
        NOconc = 0
        
    if 'NO3' in comp_namelist:
        NO3conc = y[comp_namelist.index('NO3')]
    else:
        NO3conc = 0
        
    if 'HO2' in comp_namelist:
        HO2conc = y[comp_namelist.index('HO2')]
    else:
        HO2conc = 0
        
    if 'H2O' in comp_namelist:
        H2Oconc = y[comp_namelist.index('H2O')]
    else:
        H2Oconc = 0
    return [O2conc, H2Oconc, NOconc, HO2conc, NO3conc]