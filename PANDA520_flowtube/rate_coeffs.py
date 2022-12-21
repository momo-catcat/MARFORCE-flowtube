'''module for calculating reaction rate coefficients (automatically generated)'''
# module to hold expressions for calculating rate coefficients # 
# created at 2022-12-21 17:01:47.304989

import numpy

def evaluate_rates(RO2, TEMP, lightm, M, N2, H2O,  O2, NO, HO2, NO3,p):

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
		KMT06=1+(1.40e-21*numpy.exp(2200/TEMP)*H2O) 
		w1a=4.687e-10-1.3855e-5*numpy.exp(-0.75*p/1.62265)+5.51868e-10*numpy.exp(-0.75*p/199.328) 
		w2a=-0.00331-0.00514*numpy.exp(-0.75*p/325.68711)-0.00444*numpy.exp(-0.75*p/40.81609) 
		w1b=1.1659e-9-7.79644e-10*numpy.exp(-0.75*p/22.09281)+1.03779e-9*numpy.exp(-0.75*p/568.15381) 
		w2b=-0.00813-0.00382*numpy.exp(-0.75*p/45.57591)-0.00643*numpy.exp(-0.75*p/417.95061) 

	except:
		erf = 1 # flag error
		err_mess = 'Error: reaction rates failed to be calculated, please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file' # error message

	# estimate and append photolysis rates
	J = numpy.zeros(52) 
	J[1] = 1e-5 

	rate_values = numpy.zeros((15))
	
	# reac_coef has been formatted so that python can recognize it
	# gas-phase reactions
	rate_values[0] = 0.11 * 5.4e-11 * numpy.exp(180 / TEMP)
	rate_values[1] = 0.38 * 5.4e-11 * numpy.exp(180 / TEMP)
	rate_values[2] = 0.45 * 5.4e-11 * numpy.exp(180 / TEMP)
	rate_values[3] = 2.1e-10
	rate_values[4] = w1a * numpy.exp(w2a * TEMP)
	rate_values[5] = w1b * numpy.exp(w2b * TEMP)
	rate_values[6] = 1.0e-10
	rate_values[7] = 1.6e-11 * numpy.exp(440 / TEMP)
	rate_values[8] = 2.0e-13
	rate_values[9] = 1.47e-11 * numpy.exp(- 1090 / TEMP)
	rate_values[10] = 1.4e-11 * numpy.exp(540 / TEMP)
	rate_values[11] = 2 * 6.9e-31 * (TEMP / 300) ** -0.8 * p / 1.3806488e-23 / TEMP / 1e6
	rate_values[12] = 6.2e-14 * (TEMP / 298) ** 2.6 * numpy.exp(945 / TEMP)
	rate_values[13] = 4.8e-11 * numpy.exp(250 / TEMP)
	rate_values[14] = 2.20e-13 * KMT06 * numpy.exp(600 / TEMP) + 1.90e-33 * M * KMT06 * numpy.exp(980 / TEMP)
	
	return rate_values
