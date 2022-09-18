# this script is used to calculate particle penetration rate according to
# Gormley & Kennedy, 1948, Diffusion from a stream flowing through a cylindrical tube
import numpy as np
from Vapor_calc import Air_density as Air_density
from Vapor_calc import Dp as Dp
import pandas as pd
from molmass import Formula
from thermo.chemical import Chemical


# Reynolds number of fluid
def cal_Re(dens, v, d, vis):
    """
    dens: density of air (kg m-3)
    v: velocity of air flow (m s-1)
    d: diameter of the tube (m)
    vis: viscosity of the air
    """
    return dens * v * d / vis

# Dynamic viscosity (rather than kinetic viscosity which is divided by density)
# of air using Sutherland equation, Eq. 2-8, P18, Aerosol Measurement, Third Edition
def cal_vis(T ):
    return 18.203 * 10e-6 * (293.15 + 110.4) / (T + 110.4) * (T / 293.15) ** 1.5  # Pa*s


# Mean free path of air, Eq. 2-10, P19, Aerosol Measurement, Third Edition
def cal_mfp(T, P):
    return 66.5 * 10e-9 * (101325 / P) * (T / 293.15) * (1 + 110.4 / 293.15) / (1 + 110.4 / T)  # in m


# Knudsen number
def cal_Kn(dp, mfp):
    return 2 * mfp / dp


# Cunningham slip correction factor, equation and parameters by Allen & Raabe 1985, for solid particles
# Can also be found in Eq. 2-14, P20,  Aerosol Measurement, Third Edition
def cal_Cs(Kn):
    alpha = 1.142
    beta = 0.558
    gamma = 0.999
    return 1 + Kn * (alpha + beta * np.exp(- gamma / Kn))

# Diffusion coefficient
def cal_diffusivity(dp, Cs, vis, T ):
    kB = 1.38064852e-23  # m2 kg s-2 K-1
    return kB * T * Cs / (3 * np.pi * vis * dp)  # in m2/s


# Diffusional losses: particle penetration efficiency through a Hagen-Poiseuille
# flow in a cylindrical tubing, Gormley-Kennedy formula (1949)
def cal_pene_rate(P, T, MM, rho, L, Q, Dia):
    """
    P: pressure
    T: temperature (K)
    MM: molecular weight of the molecule (Th)
    Rho: sensity of the bulk material (kg m-3) or density ???
    L: length of the inlet in m
    Q: flow rate in m3 s-1
    Dia: diameter of the tube m-1 ???
    """

    dens = Air_density(P, T)

    v = Q / (np.pi * Dia / 4)  # velocity of the gas

    vis = cal_vis(T)

    Re = cal_Re(dens, v, Dia, vis)  # Reynolds number

    dp = Dp(MM, rho) + 0.3e-9  # correction for mass diameter to mobility diameter

    mfp = cal_mfp(T, P)  # mean free path

    Kn = cal_Kn(dp, mfp)

    Cs = cal_Cs(Kn)

    D = cal_diffusivity(dp, Cs, vis, T)

    miu = np.pi * D * L / Q

    if miu > 0.02:
        pene = 0.819 * np.exp(-3.66 * miu) + 0.0975 * np.exp(-22.3 * miu) + 0.0325 * np.exp(-57 * miu) \
               + 0.0154 * np.exp(-107.6 * miu)
    else:
        pene = 1.0 - 2.56 * miu ** (2 / 3) + 1.2 * miu + 0.1767 * miu ** (4 / 3)

    export = pd.DataFrame({'Pene': [pene], 'D': [D]})

    return(export)


def cal_diffu(A,p,T):
    mass = Formula(A).mass
    
    tol = Chemical(A)
    
    tol.calculate(T=T, P=p)
    
    vis = cal_vis(T)

    dp = Dp(mass, tol.rho) + 0.3e-9  # correction for mass diameter to mobility diameter

    mfp = cal_mfp(T, p)  # mean free path

    Kn = cal_Kn(dp, mfp)

    Cs = cal_Cs(Kn)
    
    d = cal_diffusivity(dp, Cs, vis, T)
    return d*100**2