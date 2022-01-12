'''module for calculating reaction rate coefficients (automatically generated)'''
# module to hold expressions for calculating rate coefficients # 
# created at 2022-01-12 11:40:01.564803

import numpy

def evaluate_rates(RO2, H2O, TEMP, lightm, M, N2, O2, NO, HO2, NO3,p):

	# inputs: ------------------------------------------------------------------
	# RO2 - names of components included in peroxy radical list
	# M - third body concentration (molecules/cc (air))
	# N2 - nitrogen concentration (molecules/cc (air))
	# O2 - oxygen concentration (molecules/cc (air))
	# H2O, TEMP: given by the user
	# lightm: given by the user and is 0 for lights off and 1 for on
	# reaction rate coefficients and their names parsed in eqn_parser.py 
	# Jlen - number of photolysis reactions
	# tf - sunlight transmission factor
	# NO - NO concentration (# molecules/cm3 (air))
	# HO2 - HO2 concentration (# molecules/cm3 (air))
	# NO3 - NO3 concentration (# molecules/cm3 (air))
	# tf_UVC - transmission factor for 254 nm wavelength light (0-1) 
	# ------------------------------------------------------------------------

	erf = 0; err_mess = '' # begin assuming no errors
	# calculate generic reaction rate coefficients given by chemical scheme
	try:
		K120=2.5e-31*M*(TEMP/300)**-2.6 
		K12I=2.0e-12 
		KR12=K120/K12I 
		FC12=0.53 
		NC12=0.75-1.27*(numpy.log10(FC12)) 
		F12=10**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**2)) 
		KMT12=(K120*K12I*F12)/(K120+K12I) 
		KMT06=1+(1.40e-21*numpy.exp(2200/TEMP)*H2O) 

	except:
		erf = 1 # flag error
		err_mess = 'Error: reaction rates failed to be calculated, please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file' # error message

	# estimate and append photolysis rates
	J = numpy.zeros(52) 
	J[1] = 1e-5 

	rate_values = numpy.zeros((5))
	
	# reac_coef has been formatted so that python can recognize it
	# gas-phase reactions
	rate_values[0] = 1.32e-12 * (TEMP / 300) ** -0.7
	rate_values[1] = 4.8e-11*numpy.exp(250/TEMP)
	rate_values[2] = 2*6.9e-31 * (TEMP / 300) ** -0.8 * p / 1.3806488e-23 / TEMP / 1e6
	rate_values[3] = 1.3e-12*numpy.exp(-330/TEMP)
	rate_values[4] = 3.9e-41 * numpy.exp(6830.6 / TEMP)
	
	return(rate_values, erf, err_mess)
