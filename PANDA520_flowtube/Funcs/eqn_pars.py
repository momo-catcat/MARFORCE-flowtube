########################################################################
#								       #
# Copyright (C) 2018-2024					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''parses the input files to automatically create the solver file'''
# input files are interpreted and used to create the necessary
# arrays and python files to solve problem

import numpy as np
import Funcs.sch_interr as sch_interr
import Funcs.xml_interr as xml_interr
import Funcs.eqn_interr as eqn_interr
import Funcs.photo_num as photo_num 
import Funcs.write_dydt_rec as write_dydt_rec
import Funcs.write_ode_solv as write_ode_solv
import Funcs.write_rate_file as write_rate_file
import Funcs.write_hyst_eq as write_hyst_eq
import Funcs.jac_setup as jac_setup
# import aq_mat_prep

# define function to extract the chemical mechanism
def extr_mech(int_tol, num_sb, modelparams):
	
	# inputs: ----------------------------------------------------
	# modelparams.sch_name - file name of chemical scheme
	# modelparams.chem_sch_mrk - markers to identify different sections of 
	# 	the chemical scheme
	# modelparams.xml_name - name of xml file
	# modelparams.photo_path - path to file containing absorption 
	# 	cross-sections and quantum yields
	# modelparams.con_infl_nam - chemical scheme names of components with 
	# 		continuous influx
	# int_tol - integration tolerances
	# modelparams.wall_on - marker for whether to include wall partitioning
	# num_sb - number of size bins (including any wall)
	# modelparams.const_comp - chemical scheme name of components with 
	#	constant concentration
	# drh_str - string from user inputs describing 
	#	deliquescence RH (fraction 0-1) as function of temperature (K)
	# erh_str - string from user inputs describing 
	#	efflorescence RH (fraction 0-1) as function of temperature (K)
	# modelparams.dil_fac - fraction of chamber air extracted/s
	# modelparams.sav_nam - name of folder to save results to
	# modelparams.pcont - flag for whether seed particle injection is 
	#	instantaneous (0) or continuous (1)
	# modelparams - reference to PyCHAM program
	# ------------------------------------------------------------

	# # getting observation details -----------------------
	# # note this must come before the call to eqn_interr
	# # so that any seed components not included in the
	# # chemical scheme are identified	
	# # if observation file provided for constraint, then
	# # get observations
	# if hasattr(modelparams, 'obs_file') and modelparams.obs_file != []:

	# 	from obs_file_open import obs_file_open
	# 	modelparams = obs_file_open(modelparams)
	# # if no observation file, then remove any stored observed values
	# else:
	# 	if hasattr(modelparams, 'obs_comp'):
	# 		delattr(modelparams, 'obs_comp')

	# --------------------------------------------------

	# starting error flag and message (assumes no errors)
	erf = 0
	err_mess = ''
	
	f_open_eqn = open(modelparams.sch_name, mode='r') # open the chemical scheme file
	# read the file and store everything into a list
	total_list_eqn = f_open_eqn.readlines()
	f_open_eqn.close() # close file
	
	# interrogate scheme to list equations
	[rrc, rrc_name, RO2_names, modelparams] = sch_interr.sch_interr(total_list_eqn, modelparams)

	# interrogate xml to list all component names and SMILES
	[err_mess_new, modelparams] = xml_interr.xml_interr(modelparams)
	# get equation information for chemical reactions
	[comp_list, Pybel_objects, comp_num, erf, err_mess, modelparams] = eqn_interr.eqn_interr(
		num_sb, erf, err_mess, modelparams)    
	
	modelparams.RO2_names = RO2_names

	if (erf != 0): # exit and display error, if it exists
		return([], [], [], [], erf, err_mess, modelparams)                               
	                    	
	# # prepare aqueous-phase and surface (e.g. wall) reaction matrices for applying 
	# # to reaction rate calculation
	# # if aqueous-phase or surface (e.g. wall) reactions present
	# if (modelparams.eqn_num[1] > 0 or modelparams.eqn_num[2] > 0):
	# 	[] = aq_mat_prep.aq_mat_prep(num_sb, comp_num, modelparams)                                                  
	

	# # if particle-phase equations are provided by particles 
	# # not turned on then raise an error
	# if (modelparams.eqn_num[1] > 0 and num_sb-modelparams.wall_on == 0):
	# 	erf = 1 # raise error
	# 	err_mess = str('Error: ' + str(modelparams.eqn_num[1]) + ' particle-phase ' +
	# 	'reactions were registered (from the chemical scheme input file), but ' +
	# 	'no particle size bins have been invoked (number_size_bins variable in ' +
	# 	'the model variables input file). Please ensure consistency. (message ' +
	# 	'generated by eqn_pars.py module)')
	
	# get index of components with continuous influx/concentration -----------
	# empty array for storing index of components with constant influx
	# modelparams.con_infl_indx = np.zeros((len(modelparams.con_infl_nam)))
	modelparams.con_C_indx = np.zeros(len(modelparams.const_comp)).astype('int')
	delete_row_list = [] # prepare for removing rows of unrecognised components
	
	# icon = 0 # count on constant influxes

	# for i in range (len(modelparams.con_infl_nam)):
		
	# 	if (modelparams.con_infl_nam[icon] == 'H2O'): # if water influxed
	# 		# if water not included explicitly in chemical schemes 
	# 		# (note it is accounted for later in init_conc)
	# 		if (modelparams.H2O_in_cs == 2):
	# 			modelparams.con_infl_indx[icon] = int(comp_num)
	# 			icon += 1 # count on constant influxes
	# 			continue
	# 		else: # if water in chemical scheme
	# 			modelparams.con_infl_indx[icon] = modelparams.comp_namelist.index(
	# 			modelparams.con_infl_nam[icon])
	# 			icon += 1 # count on constant influxes
	# 			continue

	# 	# if we want to remove continous influxes 
	# 	# not present in the chemical scheme
	# 	if (modelparams.remove_influx_not_in_scheme == 1):

	# 		try:
	# 			# index of where components with continuous
	# 			# influx occur in list of components
	# 			modelparams.con_infl_indx[icon] = modelparams.comp_namelist.index(
	# 			modelparams.con_infl_nam[icon])
	# 		except:
	# 			# remove names of unrecognised components
	# 			modelparams.con_infl_nam = np.delete(modelparams.con_infl_nam, (icon), axis=0)
	# 			# remove emissions of unrecognised components
	# 			modelparams.con_infl_C = np.delete(modelparams.con_infl_C, (icon), axis = 0)
	# 			# remove empty indices of unrecognised components
	# 			modelparams.con_infl_indx = np.delete(
	# 			modelparams.con_infl_indx, (icon), axis = 0)
	# 			icon -= 1 # count on constant influxes
	# 	else:

	# 		try:
	# 			# index of where components with continuous 
	# 			# influx occur in list of components
	# 			modelparams.con_infl_indx[i] = modelparams.comp_namelist.index(
	# 			modelparams.con_infl_nam[i])
	# 		except:
	# 			erf = 1 # raise error
				
	# 			err_mess = str('Error: continuous influx component with name ' + str(modelparams.con_infl_nam[i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')
	
	# 	icon += 1 # count on continuous influxes

	# if len(modelparams.con_infl_nam) > 0:

	# 	# get ascending indices of components with continuous influx
	# 	si = modelparams.con_infl_indx.argsort()
	# 	# order continuous influx indices, influx 
	# 	# rates and names ascending by component index
	# 	modelparams.con_infl_indx = (modelparams.con_infl_indx[si]).astype('int')
	# 	modelparams.con_infl_C = (modelparams.con_infl_C[si])
	# 	modelparams.con_infl_nam = (modelparams.con_infl_nam[si])

	# components with constant concentration
	for i in range(len(modelparams.const_comp)):
		
		try:
			# index of where constant concentration components occur in list 
			# of components
			modelparams.con_C_indx[i] = modelparams.comp_namelist.index(
				modelparams.const_comp[i])
				
		except:
			# if a component doesn't appear at a given time then provide
			# a marker for no component
			if (modelparams.const_comp[i] == ''):
				modelparams.con_C_indx[i] = -1e6
				continue	
				
				# if water then we know it will be the next 
				# component to
				# be appended to the component list
			if (modelparams.const_comp[i] == 'H2O'):
				modelparams.con_C_indx[i] = len(
					modelparams.comp_namelist)
			else: # if not water
				erf = 1 # raise error
				err_mess = str('''Error: constant 
				concentration 
				component with name ''' + 
				str(modelparams.const_comp[i]) + ''' 
				has not been identified in the 
				chemical scheme, 
				please check it is present and the 
				chemical scheme markers are correct''')	

	# -------------------------------------------------------------
	# check if water in continuous influx components 
	# if ('H2O' in modelparams.con_infl_nam):

	# 	modelparams.H2Oin = 1 # flag for water influx

	# 	# index of water in continuous influx
	# 	wat_indx = modelparams.con_infl_nam.tolist().index('H2O')

	# 	# get influx rate of water
	# 	modelparams.con_infl_H2O = modelparams.con_infl_C[wat_indx, :].reshape(1, -1)

	# 	# do not allow continuous influx of water in the standard ode 
	# 	# solver, instead deal with it inside the water ode-solver
	# 	modelparams.con_infl_C = np.delete(modelparams.con_infl_C, wat_indx, axis=0)
	# 	modelparams.con_infl_indx = np.delete(modelparams.con_infl_indx, wat_indx, axis=0)
	# 	modelparams.con_infl_nam = np.delete(modelparams.con_infl_nam, wat_indx, axis=0)
	# else:
	# 	modelparams.H2Oin = 0 # flag for no water influx
	
	# ensure integer
	# modelparams.con_infl_indx = modelparams.con_infl_indx.astype('int')
	
	[rowvals, colptrs, modelparams] = jac_setup.jac_setup(comp_num, 0, 
		0, modelparams)
	
	# call function to generate ordinary differential equation (ODE)
	# solver module, add two to comp_num to account for water 
	# and core component
	write_ode_solv.ode_gen(int_tol, rowvals, comp_num+modelparams.H2O_in_cs, 
		0, 0, modelparams)
	# print(modelparams.reac_coef_g)
	# call function to generate reaction rate calculation module
	write_rate_file.write_rate_file(rrc, rrc_name, 0, modelparams)

	# call function to generate module that tracks change tendencies
	# of certain components
	write_dydt_rec.write_dydt_rec(modelparams)
	
	# write the module for estimating deliquescence and efflorescence 
	# relative humidities as a function of temperature
	write_hyst_eq.write_hyst_eq( modelparams)

	# get number of photolysis equations
	Jlen = photo_num.photo_num(modelparams.photo_path)

	# in case equation parsing to be skipped in further simulations
	modelparams.rowvals = rowvals; modelparams.colptrs = colptrs; modelparams.comp_num = comp_num;
	modelparams.rel_SMILES = comp_list; modelparams.Pybel_objects = Pybel_objects; modelparams.Jlen = Jlen
    
	return(rrc, rrc_name, rowvals, colptrs, comp_num, Jlen, erf, err_mess, modelparams)