## This function is used to calculate vapor properties
import numpy as np

def Air_density(P, T):
    """
    p: pressure (pa)
    T: temperature (K)
    """
    R_speci = 287.058 # specific gas constant for dry air J / (kg·K)
    return P / R_speci / T #density in kg m-3


def Dp(MM, rho):
    # This function is used to calcualte diameter of a specific molecule with a fixed density assuming a hard sphere
    # The mobility size would be the hard-sphere size + 0.3 nm
    # MM : molecular mass (Th)
    # rho : density of the particle (kg m-3)
    Na = 6.022e23
    Mass = MM / Na / 1000  # kg

    Diameter = (6 * Mass / np.pi / rho) ** (1 / 3);  # m
    return (Diameter)

