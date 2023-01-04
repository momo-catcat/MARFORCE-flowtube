def add_diff_const(p,T):
    # diffusion constants that you want to input, [cm**2/s]
    # DOH-air = 165 ± 20 Torr cm2 s-1, DHO2-He = 430 ± 30 Torr cm2 s-1，DO3-He = 410 ± 25 Torr cm2 s-1 at 296 K.
    # Source OH, HO2, and Ozone Gaseous Diffusion Coefficients
    dOH = 165 / (0.00750062 * p)  # convert to pa
    T0 = 296
    dOH = 101325 / p * dOH * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

    # Source A Measurement of the Diffusion Coefficient of Hydrogen Peroxide Vapor into Air
    dHO2 = 111 / (0.00750062 * p)  # convert to pa
    T0 = 296
    dHO2 = 101325 / p * dHO2 * ((T ** (3 / 2)) / (T0 ** (3 / 2)))

    # https://www.engineeringtoolbox.com/air-diffusion-coefficient-gas-mixture-temperature-d_2010.html
    dH2O = 0.242  # cm2/s 20C
    T0 = 273.15 + 20
    dH2O = 101325 / p * dH2O * ((T ** (3 / 2)) / (T0 ** (3 / 2)))
    return dOH, dH2O, dHO2