'''automatically produces a module for calculating reacion rate coefficients'''
# #########################################################################################
# Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk
# # All Rights Reserved.
# This file is part of PyCHAM, here, we have small modification
# # PyCHAM is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
# # PyCHAM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more details.
# # You should have received a copy of the GNU General Public License along with PyCHAM.
# If not, see <http://www.gnu.org/licenses/>.
# #########################################################################################
# function to generate a module for calculation of reaction rate coefficients

import datetime


def write_rate_file(reac_coef_g, p, rrc, rrc_name, testf):  # define function

    # inputs: ----------------------------------------------------------------------------
    # reac_coef_g - gas-phase reaction rate coefficient expression from the equation file
    # reac_coef_aq - aqueous-phase reaction rate coefficient expression from the equation file
    # rrc - expression for generic reaction rate coefficients
    # rrc_name - name given to generic reaction rate coefficients
    # testf - flag for mode: 0 in gas-phase equation mode, 2 for test mode, 3 for
    #			aqueous-phase equation mode
    # ------------------------------------------------------------------------------------

    # open/create relevant file to write module to
    if (testf == 0):
        f = open('rate_coeffs.py', mode='w')
    if (testf == 3):
        f = open('rate_coeffs_aq.py', mode='w')
    if (testf == 2):
        f = open('rate_coeffs.py', mode='w')

    f.write('\'\'\'module for calculating reaction rate coefficients (automatically generated)\'\'\'\n')
    f.write(
        '# module to hold expressions for calculating rate coefficients # \n')  # python will convert \n to os.linesep
    f.write('# created at %s\n' % (datetime.datetime.now()))
    f.write('\n')
    f.write('import numpy\n')
    # 	f.write('import photolysisRates\n')
    f.write('\n')

    # following part is the function (there should be an indent at the start of each line)
    # suggest using one tab
    f.write('def evaluate_rates(RO2, TEMP, lightm, M, N2, H2O,  O2, NO, HO2, NO3,p):\n')
    f.write('\n')
    f.write('	# inputs: ------------------------------------------------------------------\n')
    f.write('	# RO2 - names of components included in peroxy radical list\n')
    f.write('	# M - third body concentration (molecules/cc (air))\n')
    f.write('	# N2 - nitrogen concentration (molecules/cc (air))\n')
    f.write('	# O2 - oxygen concentration (molecules/cc (air))\n')
    f.write('	# H2O, TEMP: given by the user\n')
    f.write('	# lightm: given by the user and is 0 for lights off and 1 for on\n')
    f.write('	# reaction rate coefficients and their names parsed in eqn_parser.py \n')
    f.write('	# Jlen - number of photolysis reactions\n')
    f.write('	# tf - sunlight transmission factor\n')
    f.write('	# NO - NO concentration (# molecules/cm3 (air))\n')
    f.write('	# HO2 - HO2 concentration (# molecules/cm3 (air))\n')
    f.write('	# NO3 - NO3 concentration (# molecules/cm3 (air))\n')
    f.write('	# tf_UVC - transmission factor for 254 nm wavelength light (0-1) \n')
    f.write('	# ------------------------------------------------------------------------\n')
    f.write('\n')
    f.write('	erf = 0; err_mess = \'\' # begin assuming no errors')
    f.write('\n')
    f.write('	# calculate generic reaction rate coefficients given by chemical scheme\n')
    if rrc:
        f.write('	try:\n')
        # code to calculate rate coefficients given by chemical scheme file
        for line in rrc:
            f.write('		%s \n' % line)
        f.write('\n')
        f.write('	except:\n')
        f.write('		erf = 1 # flag error\n')
        f.write(
            '		err_mess = \'Error: reaction rates failed to be calculated, please check chemical scheme and associated chemical scheme markers, which are stated in the model variables input file\' # error message\n')
    f.write('\n')
    f.write('	# estimate and append photolysis rates\n')
    f.write('	J = numpy.zeros(52) \n')
    f.write('	J[1] = 1e-5 \n')
    f.write('\n')
    # 	f.write('	J[41] = 0\n')
    # 	f.write('	J[15] = 0\n')
    # 	f.write('	J[11] = 0\n')
    # 	f.write('	J[12] = 0\n')
    # 	f.write('	J[51] = 0\n')
    # in case investigating reaction rate coefficients
    # f.write('	print(KMT09, K90, F9, NC9, KR9, K9I, FC9)\n')
    # in case of indexing error when prescribing photolysis rates
    # f.write('	J[1::] = J[0:-1]; J[0] = 0.\n')

    # calculate the rate coefficient for each equation
    f.write('	rate_values = numpy.zeros((%i))\n' % (len(reac_coef_g)))
    # BE NOTIFIED!!!: before writing the script, 'reac_coef' must be converted to
    # python-compatible format
    f.write('	\n')
    f.write('	# reac_coef has been formatted so that python can recognize it\n')
    f.write('	# gas-phase reactions\n')
    for eqn_key in range(len(reac_coef_g)):
        f.write('	rate_values[%s] = %s\n' % (eqn_key, reac_coef_g[eqn_key]))
    f.write('	\n')
    # 	f.write('	# aqueous-phase reactions\n')
    # 	for eqn_key_aq in range (1, len(reac_coef_aq)+1):
    # 		f.write('	rate_values[%s] = %s\n' %(eqn_key+eqn_key_aq, reac_coef_aq[eqn_key_aq-1]))
    # 	f.write('	\n')
    f.write('	return rate_values\n')
    f.close()

    return ()
