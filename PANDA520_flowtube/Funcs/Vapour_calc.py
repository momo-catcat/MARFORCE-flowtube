import numpy as np


def H2O_conc(Ts, RH):
    '''
    The input is H2O_conc(Ts, RH)
    Ts in K
    RH in '0.1' = 10% form
    The output H2O_conc is in absolute concentration
    '''
    import numpy as np
    import pandas as pd
    # calculate saturation vapor pressure
    T = Ts
    Tc = 647.096  # K
    Pc = 220640  # hPa
    C1 = -7.85951783
    C2 = 1.84408259
    C3 = -11.7866497
    C4 = 22.6807411
    C5 = -15.9618719
    C6 = 1.80122502

    nu = 1 - T / Tc;
    sat_vapor_p_bubbler = Pc * np.exp(Tc / T * (C1 * nu + C2 * nu ** 1.5 + C3 * nu ** 3 + \
                                                C4 * nu ** 3.5 + C5 * nu ** 4 + C6 * nu ** 7.5)) * 100

    # ice
    A = 6.116441
    m = 7.591386
    Tn = 240.7263

    Pw = RH * sat_vapor_p_bubbler
    DPs = Tn / (m / (np.log10(Pw / 100 / A)) - 1)
    H2Oconcs = sat_vapor_p_bubbler * RH / 1.3806488e-23 / T / 1e6  # #/cm3

    # define value
    Water_datalib = {'DPs': [DPs], 'H2O_conc': [H2Oconcs], 'SatP': sat_vapor_p_bubbler}
    Water = pd.DataFrame(Water_datalib, columns=['DPs', 'H2O_conc', 'SatP'])
    return (Water)


def Air_density(P, T):
    """
    p: pressure (pa)
    T: temperature (K)
    """
    R_speci = 287.058  # specific gas constant for dry air J / (kg·K)
    return P / R_speci / T  # density in kg m-3


def Dp(MM, rho):
    # This function is used to calcualte diameter of a specific molecule with a fixed density assuming a hard sphere
    # The mobility size would be the hard-sphere size + 0.3 nm
    # MM : molecular mass (Th)
    # rho : density of the particle (kg m-3)
    Na = 6.022e23
    Mass = MM / Na / 1000  # kg

    Diameter = (6 * Mass / np.pi / rho) ** (1 / 3);  # m
    return (Diameter)
