'''module to calculate gas-phase concentration of water'''
# based on input relative humidity and temperature the
# concentration of water in the gas-phase is estimated

import numpy as np
import scipy.constants as si
def water_calc(TEMP, RH):
    
	# saturation vapour pressure (Pa) (TEMP in K.)
	Psat_water = (np.exp((-0.58002206E4/TEMP)+0.13914993E1-(0.48640239E-1*TEMP) 
		+ (0.41764768E-4*(TEMP**2e0))-(0.14452093E-7*(TEMP**3e0))+ 
		(0.65459673e1*np.log(TEMP))))
	NA = si.Avogadro # Avogadro's number (molecules/mol)

	# convert saturation vapour pressures from Pa to molecules/cc (air) using ideal
	# gas law, R has units cc.Pa/K.mol
	H2Oc = (RH*Psat_water)*(NA/(8.3144598e6*TEMP))
	
	# convert Psat_water to log10(atm) from Pa
	Psat_water = np.log10(Psat_water/101325.0)
	H2O_mw = 18.0 # state molecular weight of water (g/mol)
	
	return(H2Oc, Psat_water, H2O_mw)
