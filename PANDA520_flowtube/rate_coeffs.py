'''module for calculating reaction rate coefficients (automatically generated)'''
# module to hold expressions for calculating rate coefficients # 
# created at 2022-01-31 11:09:43.362268

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

	except:
		erf = 1 # flag error
		err_mess = 'Error: reaction rates failed to be calculated, please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file' # error message

	# estimate and append photolysis rates
	J = numpy.zeros(52) 
	J[1] = 1e-5 

	rate_values = numpy.zeros((3))
	
	# reac_coef has been formatted so that python can recognize it
	# gas-phase reactions
	rate_values[0] = 1.8e-10
	rate_values[1] = 2 * 6.9e-31 * (TEMP / 300) ** -0.8 * p / 1.3806488e-23 / TEMP / 1e6
	rate_values[2] = 4.8e-11 * numpy.exp(250 / TEMP)
	
	return(rate_values, erf, err_mess)
