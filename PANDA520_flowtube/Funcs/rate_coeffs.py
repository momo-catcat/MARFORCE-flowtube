##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2024 Simon O'Meara : simon.omeara@manchester.ac.uk               #
#                                                                                        #
#    All Rights Reserved.                                                                #
#    This file is part of PyCHAM                                                         #
#                                                                                        #
#    PyCHAM is free software: you can redistribute it and/or modify it under             #
#    the terms of the GNU General Public License as published by the Free Software       #
#    Foundation, either version 3 of the License, or (at your option) any later          #
#    version.                                                                            #
#                                                                                        #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #
#    details.                                                                            #
#                                                                                        #
#    You should have received a copy of the GNU General Public License along with        #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #
#                                                                                        #
##########################################################################################
'''module for calculating reaction rate coefficients (automatically generated)'''
# module to hold expressions for calculating rate coefficients # 
# created at 2026-07-06 14:19:35.118903

import numpy
from Funcs.Wennberg_rec_funcs import TROE, TUN, ALK, NIT, EPO, ISO1, ISO2, KCO
from Funcs.itprod import itProd
from Funcs.kclust import kClust, kDimer, kTrimer

def evaluate_rates(RO2, H2O, TEMP, time, M, N2, O2, Jlen, NO, HO2, NO3, sumt, p, modelparams):

	# inputs: ------------------------------------------------------------------
	# RO2 - total concentration of alkyl peroxy radicals (# molecules/cm3) 
	# M - third body concentration (# molecules/cm3 (air))
	# N2 - nitrogen concentration (# molecules/cm3 (air))
	# O2 - oxygen concentration (# molecules/cm3 (air))
	# H2O, TEMP: given by the user
	# modelparams.light_stat_now: given by the user and is 0 for lights off and >1 for on
	# reaction rate coefficients and their names parsed in eqn_parser.py 
	# Jlen - number of photolysis reactions
	# modelparams.tf - sunlight transmission factor
	# NO - NO concentration (# molecules/cm3 (air))
	# HO2 - HO2 concentration (# molecules/cm3 (air))
	# NO3 - NO3 concentration (# molecules/cm3 (air))
	# modelparams.tf_UVC - transmission factor for 254 nm wavelength light (0-1) 
	# ------------------------------------------------------------------------

	erf = 0; err_mess = '' # begin assuming no errors
	SUN = getattr(modelparams, "sun_a", 0)
	print('SUN:', SUN)

	# calculate any generic reaction rate coefficients given by chemical scheme

	try:
		gprn=0
		gprn += 1 # keep count on reaction number
		KMT06=1+(1.40e-21*numpy.exp(2200/TEMP)*H2O); 
		gprn += 1 # keep count on reaction number
		K120=2.5e-31*M*(TEMP/300)**(-2.6); 
		gprn += 1 # keep count on reaction number
		K12I=2.0e-12; 
		gprn += 1 # keep count on reaction number
		KR12=K120/K12I; 
		gprn += 1 # keep count on reaction number
		FC12=0.53; 
		gprn += 1 # keep count on reaction number
		NC12=0.75-1.27*(numpy.log10(FC12)); 
		gprn += 1 # keep count on reaction number
		F12=10**(numpy.log10(FC12)/(1.0+(numpy.log10(KR12)/NC12)**(2))); 
		gprn += 1 # keep count on reaction number
		KMT12=(K120*K12I*F12)/(K120+K12I); 
		gprn += 1 # keep count on reaction number
		KMT05=1.44e-13*(1+(M/4.2e+19)); 

	except:
		erf = 1 # flag error
		err_mess = str('Error: generic reaction rates failed to be calculated inside rate_coeffs.py at number ' + str(gprn) + ', please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file') # error message
		return([], erf, err_mess)

	if (modelparams.light_stat_now == 0):
		J = [0]*Jlen
	rate_values = numpy.zeros((15))
	
	# if reactions have been found in the chemical scheme
	# gas-phase reactions
	gprn = 0 # keep count on reaction number
	try:
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.32e-12*(TEMP/300)**-0.7' 
		rate_values[0] = 1.32e-12*(TEMP/300)**-0.7
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '4.8e-11*numpy.exp(250/TEMP)' 
		rate_values[1] = 4.8e-11*numpy.exp(250/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '2.20e-13*KMT06*numpy.exp(600/TEMP)+1.90e-33*M*KMT06*numpy.exp(980/TEMP)' 
		rate_values[2] = 2.20e-13*KMT06*numpy.exp(600/TEMP)+1.90e-33*M*KMT06*numpy.exp(980/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.9e-31*(TEMP/300)**-0.8*p/1.3806488e-23/TEMP/1e6' 
		rate_values[3] = 6.9e-31*(TEMP/300)**-0.8*p/1.3806488e-23/TEMP/1e6
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '6.2e-14*(TEMP/298)**2.6*numpy.exp(945/TEMP)' 
		rate_values[4] = 6.2e-14*(TEMP/298)**2.6*numpy.exp(945/TEMP)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '1.3e-12*numpy.exp(-330/TEMP)*O2' 
		rate_values[5] = 1.3e-12*numpy.exp(-330/TEMP)*O2
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '3.9e-41*numpy.exp(6830.6/TEMP)*H2O*H2O' 
		rate_values[6] = 3.9e-41*numpy.exp(6830.6/TEMP)*H2O*H2O
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'kDimer(modelparams)' 
		rate_values[7] = kDimer(modelparams)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'kTrimer(modelparams)' 
		rate_values[8] = kTrimer(modelparams)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'kTrimer(modelparams)' 
		rate_values[9] = kTrimer(modelparams)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'kTrimer(modelparams)' 
		rate_values[10] = kTrimer(modelparams)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'kTrimer(modelparams)' 
		rate_values[11] = kTrimer(modelparams)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'KMT05' 
		rate_values[12] = KMT05
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = 'itProd(modelparams)' 
		rate_values[13] = itProd(modelparams)
		gprn += 1 # keep count on reaction number
		# remember equation in case needed for error reporting
		rc_eq_now = '5.4e-32*(TEMP/300)**-1.8*N2' 
		rate_values[14] = 5.4e-32*(TEMP/300)**-1.8*N2
	except:
		erf = 1 # flag error
		err_mess = (str('Error: Could not calculate rate coefficient for equation number ' + str(gprn) + ' ' + rc_eq_now + ' (message from rate coeffs.py)'))
	
	# aqueous-phase reactions
	
	# surface (e.g. wall) reactions
	
	return(rate_values, erf, err_mess)
