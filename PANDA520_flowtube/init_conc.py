'''function to initiate concentrations of components'''
# based on inputs, initial concentrations and their holding arrays 
# are set

import numpy as np
import scipy.constants as si
import math
from water_calc import water_calc
import write_dydt_rec


def init_conc(num_comp, Comp0, init_conc, TEMP, PInit, 
	 dydt_trak, rindx, pindx, num_eqn, nreac, nprod, 
	comp_namelist, RO2_indx, HOMRO2_indx, rstoi, pstoi):
# def init_conc(num_comp, Comp0, init_conc, TEMP, RH, PInit, Pybel_objects,
# 	testf, pconc, dydt_trak, end_sim_time, save_step, 
# 	rindx, pindx, num_eqn, nreac, nprod, 
# 	comp_namelist, Compt, seed_name, seed_mw,
# 	core_diss, nuc_comp, comp_xmlname, comp_smil, rel_SMILES,
# 	RO2_indx, HOMRO2_indx, rstoi, pstoi)		
	# inputs:------------------------------------------------------
	
	# num_comp - number of unique components
	# Comp0 - chemical scheme names of components present at start of experiment
	# init_conc - initial concentrations of components (ppb)	
	# TEMP - temperature in chamber at start of experiment (K)
	# RH - relative humidity in chamber (dimensionless fraction 0-1)
	# PInit - initial pressure (Pa)
	# init_SMIL - SMILES of components present at start of experiment (whose 
	# concentrations are given in init_conc)
	# testf - flag for whether in normal mode (0) or testing mode (1/2)
	# pconc - initial concentration of particles (# particles/cc (air))
	# dydt_trak - chemical scheme name of components for which user wants the tendency to  
	#			change tracked
	# end_sim_time - total simulation time (s)
	# save_step - recording frequency (s)
	# rindx - indices of reactants per equation
	# pindx - indices of products per equation
	# num_eqn - number of equations
	# comp_namelist - list of names of components as presented in the chemical scheme file
	# Compt - name of components injected instantaneously after start of experiment
	# seed_name - name of core component (input by user)
	# seed_mw - molecular weight of seed material (g/mol)
	# core_diss - dissociation constant of seed material
	# nuc_comp - name of nucleating component (input by user, or defaults to 'core')
	# comp_xmlname - component names in xml file
	# comp_smil - all SMILES strings in xml file
	# rel_SMILES - only the SMILES strings of components present in the chemical scheme file
	# RO2_indx - RO2 list indices and chemical scheme indices of non-HOM-RO2 molecules
	# HOMRO2_indx - chemical scheme indices of HOM-RO2 molecules
	# rstoi - stoichiometry of reactants per equation
	# pstoi - stoichiometry of products per equation
	# -----------------------------------------------------------

	# start by assuming no error
	erf = 0
	err_mess = ''
	
# 	if testf==1: # testing mode
# 		# return dummies
# 		return(0,0,0,0,0,0,0,0)

	NA = si.Avogadro # Avogadro's number (molecules/mol)
	# empty array for storing species' concentrations, must be an array
	y = np.zeros((num_comp))
	y_mw = np.zeros((num_comp, 1)) # species' molecular weight (g/mol)
	# empty array for storing index of interesting gas-phase components
	y_indx_plot = []
	
	# convert concentrations
	# total number of molecules in 1 cc air using ideal gas law.  R has units cc.Pa/K.mol
	ntot = PInit*(NA/((si.R*1.e6)*TEMP))
	# one billionth of number of # molecules in chamber unit volume
	Cfactor = ntot*1.e-9 # ppb to # molecules/cm3 conversion factor

	# prepare dictionary for tracking tendency to change of user-specified components
	dydt_vst = {}

	# insert initial concentrations where appropriate
	for i in range(len(Comp0)):
    		# index of where initial components occur in list of components
		try: # in case components already listed via interpretation of the chemical scheme
			y_indx = comp_namelist.index(Comp0[i])
			
		# if component not already listed via interpretation of the chemical scheme
		# then send error message
		except:
			erf = 1
			err_mess = str('Error: component called ' + str(Comp0[i]) + ', which has an initial concentration specified in the model variables input file has not been found in the chemical scheme.  Please check the scheme and associated chemical scheme markers, which are stated in the model variables input file.')
			return (0, 0, 0, 0, 0, 0, 0, 0, 
				0, 0, 0,
				0, 0, 0, 0, erf, err_mess, 0, 0, 0)
			
		y[y_indx] = init_conc[i] #*Cfactor # convert from ppb to # molecules/cm3 (air)
		
		# remember index for plotting gas-phase concentrations later
		y_indx_plot.append(y_indx)

	# number of recording steps
# 	nrec_steps = int(math.ceil(end_sim_time/save_step)+1)
	
# 	for i in range(num_comp): # loop through all components to get molar weights
# 		y_mw[i] = Pybel_objects[i].molwt # molecular weight (g/mol)
	
	# ------------------------------------------------------------------------------------
	# account for water's properties
	
	# get initial gas-phase concentration (molecules/cc (air)) and vapour pressure
	# of water (log10(atm))
# 	[C_H2O, Psat_water, H2O_mw] = water_calc(TEMP, RH[0])
	
	# holder for water index (will be used if not identified in chemical scheme)
# 	H2Oi = num_comp # index for water
# 	
# 	# check for water presence in chemical scheme via its SMILE string
# 	# count on components
# 	indx = -1
# 	for single_chem in rel_SMILES:
# 		indx += 1
# 		# ensure this is water rather than single oxygen (e.g. due to ozone photolysis 
# 		# (O is the MCM chemical scheme name for single oxygen))
# 		if (single_chem == 'O' and comp_namelist[indx] != 'O'):
# 			H2Oi = indx
# 			y[H2Oi] = C_H2O # include initial concentration of water (molecules/cm3)
# 			y_mw[H2Oi] = H2O_mw # include molar weight of water (g/mol)
	
	# if not included in chemical scheme file, then add water to end of component list
# 	if (H2Oi == num_comp):
# 	
# 		num_comp += 1 # update number of components to account for water
# 		# append empty element to y and y_mw to hold water values
# 		y = np.append(y, C_H2O)
# 		# append molar weight of water (g/mol)
# # 		y_mw = (np.append(y_mw, H2O_mw)).reshape(-1, 1)
# 		comp_namelist.append('H2O') # append water's name to component name list

	# ------------------------------------------------------------------------------------
	# account for seed properties - note that even if no seed particle, this code ensures
	# that an index is provided for core material

	# empty array for index of core component
# 	seedi = (np.zeros((len(seed_name)))).astype(int)
	
# 	comp_namelist.append('core') # append name of core to component name list
# 	corei = [num_comp] # index for core component
	# increase number of components to account for 'core' component
# 	num_comp += 1

	# append core gas-phase concentration (molecules/cc (air)) and molecular 
	# weight (g/mol) (needs to have a 1 length in second dimension for the kimt 
	# calculations)
# 	y = np.append(y, 0.)
# 	y_mw = (np.append(y_mw, seed_mw)).reshape(-1, 1)
	
	# --------------------------------------
	# get index of user-specified components for tracking their change tendencies (dydt) due to modeled
	# mechanisms
	if (len(dydt_trak) > 0):
		
		dydt_traki = [] # empty list for indices of these components
		
		for i in range (len(dydt_trak)):

			reac_index = [] # indices of reactions involving this component
			# value for whether this component is reactant or product in a reaction and stoichiometry
			reac_sign = []
			
			if dydt_trak[i] != 'RO2' and dydt_trak[i] != 'HOMRO2':
			
				# index of components in component list
				try:
					y_indx = comp_namelist.index(dydt_trak[i])
				# if component not already listed via interpretation of the chemical scheme
				# then send error message
				except:
					erf = 1
					err_mess = str('Error: component called ' + str(dydt_trak[i]) + ', which is specified to be tracked in the model variables input file has not been found in the chemical scheme.  Please check the scheme and associated chemical scheme markers, which are stated in the model variables input file.')
					return (0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0,
					0, 0, 0, 0, erf, err_mess, 0, 0, 0)
				# remember index for plotting gas-phase concentrations later
				dydt_traki.append([int(y_indx)])
			
				# search through reactions to see where this component is reactant or product
				for ri in range(num_eqn):
					if sum(rindx[ri, 0:nreac[ri]] == y_indx) > 0:
						reac_index.append(int(ri)) # append reaction index
						reac_place = np.where(rindx[ri, 0:nreac[ri]] == y_indx)[0]
						reac_sign.append(-1*rstoi[int(ri), reac_place])
					if sum(pindx[ri, 0:nprod[ri]] == y_indx) > 0:
						reac_index.append(int(ri)) # append reaction index
						reac_place = np.where(pindx[ri, 0:nprod[ri]] == y_indx)[0]
						reac_sign.append(1*pstoi[int(ri), reac_place])
			
			if dydt_trak[i] == 'RO2': # for non-HOM-RO2 (organic peroxy radicals)
			
				# remember index for plotting gas-phase concentrations later
				dydt_traki.append(list(RO2_indx[:, 1]))
				
				# loop through non-HOM-RO2 components to get their reaction indices
				for y_indx in RO2_indx[:, 1]:
					
					# search through reactions to see where this component is reactant or product
					for ri in range(num_eqn):
						if sum(rindx[ri, 0:nreac[ri]] == y_indx) > 0:
							reac_index.append(int(ri)) # append reaction index
							reac_place = np.where(rindx[ri, 0:nreac[ri]] == y_indx)[0]
							reac_sign.append(-1*rstoi[int(ri), reac_place])
							
						if sum(pindx[ri, 0:nprod[ri]] == y_indx) > 0:
							reac_index.append(int(ri)) # append reaction index
							reac_place = np.where(pindx[ri, 0:nprod[ri]] == y_indx)[0]
							reac_sign.append(1*pstoi[int(ri), reac_place])
					
					y_indx = RO2_indx[:, 1] # ready for storing below
			
			if dydt_trak[i] == 'HOMRO2': # for HOM-RO2 (highly oxidised molecule organic peroxy radicals)
			
				# remember index for plotting gas-phase concentrations later
				dydt_traki.append(list(np.squeeze(HOMRO2_indx)))
				
				# loop through non-HOM-RO2 components to get their reaction indices
				for y_indx in HOMRO2_indx[:]:
					
					# search through reactions to see where this component is reactant or product
					for ri in range(num_eqn):
						if sum(rindx[ri, 0:nreac[ri]] == y_indx) > 0:
							reac_index.append(int(ri)) # append reaction index
							reac_place = np.where(rindx[ri, 0:nreac[ri]] == y_indx)[0]
							reac_sign.append(-1*rstoi[int(ri), reac_place])
						if sum(pindx[ri, 0:nprod[ri]] == y_indx) > 0:
							reac_index.append(int(ri)) # append reaction index
							reac_place = np.where(pindx[ri, 0:nprod[ri]] == y_indx)[0]
							reac_sign.append(1*pstoi[int(ri), reac_place])
				y_indx = HOMRO2_indx[:] # ready for storing below
				
			# save reaction indices in dictionary value for this component,
			# when creating empty rec_array, add two columns onto the end for particle- and 
			# wall-partitioning, respectively
			rec_array = np.zeros((2, len(reac_index)))
			rec_array[0, :] = reac_index

			comp_indx_str = str(dydt_trak[i] + '_comp_indx')
			res_string = str(dydt_trak[i] + '_res')
			reac_string = str(dydt_trak[i] + '_reac_sign')
			
			dydt_vst[comp_indx_str] = y_indx # dictionary entry to hold index of tracked component
			dydt_vst[res_string] = rec_array # dictionary entry to hold reaction indices and results
			dydt_vst[reac_string] = reac_sign # dictionary entry to hold results


		# dictionary entry to hold component names of components to track
		dydt_vst['comp_names'] = dydt_trak
		
		# call on write_dydt_rec to generate the module that will process
		# the tendency to change during the simulation
		write_dydt_rec.write_dydt_rec()
	
	# --------------------------------------
	
# 	# if nucleating component formed of core component
# 	if (nuc_comp[0] == 'core'):
# 		nuci = num_comp-1 # index of core component
# 	else:
# 		nuci = -1 # filler
	
	# get indices of seed particle component(s)
# 	indx = 0 # count on seed component(s)
# 	for sname in seed_name:
# 		# index of core component
# 		seedi[indx] = int(comp_namelist.index(sname))
# 		indx += 1 # count on seed component(s)
# 	
# 	# get index of component with latter injections
# 	if len(Compt)>0:
# 		inj_indx = np.zeros((len(Compt)))
# 		for i in range(len(Compt)):
# 			# index of where instantaneously injected components 
# 			# occur in SMILES string
# 			inj_indx[i] = comp_namelist.index(Compt[i])
# 	else:
# 		inj_indx = np.zeros((1)) # dummy

	# ensure index arrays are integer type
# 	inj_indx = inj_indx.astype('int')
# 	corei = np.array((corei)).astype('int')
	
	# get indices of NO, HO2 and NO3 (for reaction rate calculations)
	
# 	try:
# 		NOi = comp_namelist.index('NO')
# 	except:
# 		NOi = 0 # filler
# 	try:
# 		HO2i = comp_namelist.index('HO2')
# 	except:
# 		HO2i = 0 # filler
# 	try:
# 		NO3i = comp_namelist.index('NO3')
# 	except:
# 		NO3i = 0 # filler
# 	RO2i = comp_namelist.index('RO2')
	return (y,  num_comp, Cfactor*1e9, dydt_vst,
			comp_namelist)
