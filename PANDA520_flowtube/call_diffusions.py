# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 16:13:10 2022

@author: jiali
"""
# this file is used for call diffusion for all species 
import pene_rate
from molmass import Formula
from thermo.chemical import Chemical


def cal_diffu(A,p,T):
    mass = Formula(A).mass
    
    tol = Chemical(A)
    
    tol.calculate(T=T, P=p)
    
    d = pene_rate.cal_diffu(T, p, mass, tol.rho)
    return d
            


