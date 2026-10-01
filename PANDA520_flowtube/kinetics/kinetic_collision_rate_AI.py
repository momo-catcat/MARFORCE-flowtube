"""
Kinetic collision rate calculator using gas kinetic theory.

Calculates the hard-sphere bimolecular collision rate constant:
    k_coll = pi * (r_A + r_B)^2 * sqrt(8 * kB * T / (pi * mu))

where:
    r_A, r_B  = molecular radii estimated from liquid density
    mu        = reduced mass = m_A * m_B / (m_A + m_B)
    kB        = Boltzmann constant
    T         = temperature (K)

Molecular radius is estimated from:
    r = (3 * M / (4 * pi * rho * N_A))^(1/3)

Reference: Seinfeld & Pandis, Atmospheric Chemistry and Physics, Chapter 9
           Fuller et al. (1965, 1966, 1969) for diffusion volumes
"""

import numpy as np

# Constants
kB = 1.381e-23   # Boltzmann constant, J/K
NA = 6.022e23    # Avogadro's number, mol-1


def molecular_radius(M, rho):
    """
    Estimate molecular radius from molecular weight and liquid density.

    Parameters
    ----------
    M : float
        Molecular weight (g/mol)
    rho : float
        Liquid density (g/cm3)

    Returns
    -------
    r : float
        Molecular radius (cm)
    """
    V_mol = M / (rho * NA)  # molecular volume in cm3
    r = (3 * V_mol / (4 * np.pi))**(1/3)
    return r


def collision_rate(M_A, M_B, r_A, r_B, T):
    """
    Calculate hard-sphere bimolecular collision rate constant.

    Parameters
    ----------
    M_A, M_B : float
        Molecular weights (g/mol)
    r_A, r_B : float
        Molecular radii (cm)
    T : float
        Temperature (K)

    Returns
    -------
    k_coll : float
        Collision rate constant (cm3/molec/s)
    """
    m_A = M_A * 1e-3 / NA  # kg per molecule
    m_B = M_B * 1e-3 / NA
    mu = (m_A * m_B) / (m_A + m_B)  # reduced mass (kg)

    sigma = np.pi * (r_A + r_B)**2  # collision cross section (cm2)
    v_rel = np.sqrt(8 * kB * T / (np.pi * mu)) * 100  # mean relative speed (cm/s)

    k_coll = sigma * v_rel  # cm3/molec/s
    return k_coll


if __name__ == "__main__":
    T = 293.15  # K
    rho_SA = 1.84  # g/cm3, liquid H2SO4 density

    # SA monomer: H2SO4, M = 98 g/mol
    M_mono = 98.0
    r_mono = molecular_radius(M_mono, rho_SA)

    # SA dimer: (H2SO4)2, M = 196 g/mol
    M_dimer = 196.0
    r_dimer = molecular_radius(M_dimer, rho_SA)

    # SA trimer: (H2SO4)3, M = 294 g/mol
    M_trimer = 294.0
    r_trimer = molecular_radius(M_trimer, rho_SA)

    print("=" * 65)
    print("  Kinetic collision rates for H2SO4 clustering")
    print("  T = {:.2f} K, rho(liquid H2SO4) = {:.2f} g/cm3".format(T, rho_SA))
    print("=" * 65)

    print(f"\n  Molecular radii:")
    print(f"    SA monomer: {r_mono*1e8:.2f} Angstrom")
    print(f"    SA dimer:   {r_dimer*1e8:.2f} Angstrom")
    print(f"    SA trimer:  {r_trimer*1e8:.2f} Angstrom")

    reactions = [
        ("SA + SA",         M_mono,   M_mono,   r_mono,   r_mono),
        ("SA + dimer",      M_mono,   M_dimer,  r_mono,   r_dimer),
        ("SA + trimer",     M_mono,   M_trimer, r_mono,   r_trimer),
        ("dimer + dimer",   M_dimer,  M_dimer,  r_dimer,  r_dimer),
        ("dimer + trimer",  M_dimer,  M_trimer, r_dimer,  r_trimer),
        ("trimer + trimer", M_trimer, M_trimer, r_trimer, r_trimer),
    ]

    print(f"\n  {'Reaction':<20} {'k_coll (cm3/s)':>16} {'k_coll/4e-14':>14}")
    print(f"  {'-'*52}")

    for name, ma, mb, ra, rb in reactions:
        k = collision_rate(ma, mb, ra, rb, T)
        print(f"  {name:<20} {k:>16.2e} {k/4e-14:>14.0f}x")

    print(f"\n  Note: H2SO4 has a large dipole moment (~2.7 D).")
    print(f"  Dipole-dipole attraction can enhance rates by 2-10x.")
    print(f"  kClust = 4e-14 implies sticking probability ~{4e-14/collision_rate(M_mono, M_mono, r_mono, r_mono, T):.1e}")
    print("=" * 65)
