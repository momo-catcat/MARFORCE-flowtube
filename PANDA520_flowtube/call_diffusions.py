# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 16:13:10 2022

@author: jiali
"""
# this file is used for call diffusion for all species 
from pene_rate import cal_pene_rate as cal_pene_rate
from molmass import Formula
from thermo.chemical import Chemical


def call_diff(comp_namelist,p,T,L1,Q1,R1):
    d_tot = []
    for i in range(len(comp_namelist)):
        mass = Formula(comp_namelist[i]).mass
        tol = Chemical(comp_namelist[i])
        tol.calculate(T=T, P=p)
        d = cal_pene_rate(p, T, mass, tol.rhog, L1/100, Q1/1000*60, R1/100*2)
        d_tot = [d_tot,d]

    



