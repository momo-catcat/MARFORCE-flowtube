# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 17:10:03 2022

@author: jiali
"""
# get the diffusion for all species 
from kinetics.diff_coef import diff_coef


def get_diff_and_u(comp_namelist,Diff_setname,con_C_indx,Diff_set,T,p):
    
    Diff_vals = [diff_coef(i, T, p, 'air') for i in comp_namelist]
  
    for i in Diff_setname:
        Diff_vals[comp_namelist.index(i)] = Diff_set[Diff_setname.index(i)]
    
    u = list(range(len(comp_namelist))) #% u is the index of species in C except constant compounds 
    u = [i for j, i in enumerate(u) if j not in con_C_indx]
    
    
    return(u, Diff_vals)



