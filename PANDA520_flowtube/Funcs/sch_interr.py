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
'''isolate sections of chemical scheme'''
# using the chemical scheme markers, sections of the 
# chemical scheme are separated

import re
import numpy as np
import Funcs.formatting as formatting

# define function
def sch_interr(total_list_eqn, modelparams):

	# inputs: ------------------------------------------------------------------
	# total_list_eqn - all lines from the chemical scheme file
	# modelparams - reference to PyCHAM
	# modelparams.chem_sch_mrk - markers to denote different section of the chemical scheme
	# --------------------------------------------------------------------------
	
	# preparatory part ---------------------------------------------------------
	modelparams.eqn_list = [] # empty list for gas-phase reaction equation
	modelparams.aqeqn_list = [] # empty list for particle-phase reaction equation
	modelparams.sueqn_list = [] # empty list for surface (e.g. wall) reaction equation
	RO2_names = [] # empty list for peroxy radicals
	rrc = [] # empty list for reaction rate coefficients
	rrc_name = [] # empty list for reaction rate coefficient labels
	eqn_flag = 0 # don't collate reaction equations until seen
	pr_flag = 0 # don't collate peroxy radicals until seen
	RO2_count = 0 # count on number of lines considered in peroxy radical list
	# if thirtheenth marker missing, assume this is due to old inputs
	# when only 12 markers were needed
	if (len(modelparams.chem_sch_mrk) == 12):
		modelparams.chem_sch_mrk.append('')

	# -------------------------------------------------------------------------
	
	# obtain lists for reaction rate coefficients, peroxy radicals 
	# and equation reactions using markers for separating chemical scheme elements
	for line in total_list_eqn:
		
		line1 = line.strip() # remove bounding white space
		
		# --------------------------------------------------------------------------------
		# generic reaction rate coefficients part
		# marker at end of generic reaction rate coefficient lines
		# the first \ allows python to interpret the second \ as a dash
		# to use in regex which means an escape in case the marker is a
		# regex special character
		# the $ means occurs at end of string
		end_mark = str('\\' + modelparams.chem_sch_mrk[7]+'$')
		# look out for start of generic reaction rate coefficients
		# could be generic reaction coefficient if just one = in line
	
		if (len(line1.split('=')) == 2):
			rrc_flag = 1
			# don't record if nothing preceding '=' (can occur in KPP files, e.g.
			# =IGNORE)

			if (len((line1.split('=')[0]).strip()) == 0):
				rrc_flag = 0
			# don't record if this just an IGNORE command
			if len((line1.split('=')[1]).strip()) >= 6:
				if (line1.split('=')[1]).strip()[0:6] == 'IGNORE':
					rrc_flag = 0
			if len((line1.split('=')[0]).strip()) >= 10:
				rrc_flag = 0
			# don't record if marker (if one present) for end of generic reaction rate
			# coefficient lines not present
			if (len(modelparams.chem_sch_mrk[7]) > 0):
				if re.search(end_mark, line1.strip()) == None:
					rrc_flag = 0
			
			if (rrc_flag == 1): 
				# dont consider if start of peroxy radical list
				if (line1.split('=')[0]).strip() != modelparams.chem_sch_mrk[1]:
					# don't consider if a gas-phase chemical scheme reaction
					if ((line1.split('=')[0]).strip())[0] != modelparams.chem_sch_mrk[-1]:
						# don't consider if an aqueous-phase chemical scheme reaction
						if ((line1.split('=')[0]).strip())[0] != modelparams.chem_sch_mrk[8]:
							# don't consider if a surface (e.g. wall) chemical scheme reaction
							if ((line1.split('=')[0]).strip())[0] != modelparams.chem_sch_mrk[12]:
								# remove end characters
								line2 = line1.replace(str(modelparams.chem_sch_mrk[7]), '')
								# remove all white space
								line2 = line2.replace(' ', '')
								# convert fortran-type scientific notation to python type
								line2 = formatting.SN_conversion(line2)
								# ensure rate coefficient is python readable
								line2 = formatting.convert_rate_mcm(line2)
								rrc.append(line2.strip())
								# get just name of generic reaction rate coefficient
								rrc_name.append((line2.split('=')[0]).strip())	

		# --------------------------------------------------------------------------------
		# peroxy radical part
		# start logging peroxy radicals
		RO2_start_mark = str(modelparams.chem_sch_mrk[1]+' =')
		# # Check if the starting marker for RO2 is seen
		  # Remove leading/trailing spaces

		if line1[-1] == modelparams.chem_sch_mrk[2]:
		# 	# Remove the `RO2 =` prefix if present
			if "=" in line1:
				line1 = line1.split("=", 1)[1]

		# 	# Handle continuation across multiple lines
			if ";" in line1:  # End of the RO2 list
				line1 = line1.split(";")[0]  # Remove the trailing `;`

		# 	# Split the line by the peroxy radical separator (`+`)
			species = [s.strip() for s in line1.split(modelparams.chem_sch_mrk[2])]
			for specie in species:
				# Filter out empty or invalid entries
				if specie and specie != '&':
					RO2_names.append(specie)

		elif line1[-1] == ';':
			if (re.search(r'\+', line1)) and all(c not in line1 for c in ['=', '%', '*']):
				line1 = line1.split(";")[0]  # Remove the trailing `;`
				species = [s.strip() for s in line1.split(modelparams.chem_sch_mrk[2])]
				for specie in species:
					# Filter out empty or invalid entries
					if specie and specie != '&':
						RO2_names.append(specie)
		#  if starting marker for peroxy radical list seen, flag that recording needed
		# if (re.match(RO2_start_mark, line1) != None):
			
		# 	# to double check that recording needed for peroxy radicals (in case 
		# 	# modelparams.chem_sch_mrk[1] is not unique)
			
		# 	# first check whether the RO2 list comprises just one line, as this will
		# 	# mean its end marker is present
		# 	if (len(modelparams.chem_sch_mrk[9].strip()) > 0):
		# 		# .* allows search across all elements of line, \\ ensures marker is
		# 		# recognised as string
		# 		mark = str('.*\\' + modelparams.chem_sch_mrk[9])
		# 		if (re.match(mark, line1) != None):
		# 			pr_flag = 1
			
		# 	# look for presence of marker for RO2 list continuing onto next line, which
		# 	# confirms this is the RO2 list when it covers more than one line
		# 	# .* allows search across all elements of line, \\ ensures marker is
		# 	# recognised as string
		# 	mark = str('.*\\' + modelparams.chem_sch_mrk[2])
		# 	print(re.match(mark, line1))
		# 	if (re.match(mark, line1) != None):
		# 		pr_flag = 1
			
		# 	# if line end or continuation marker not supplied then assume the RO2 start
		# 	# marker is unique
		# 	if (len(modelparams.chem_sch_mrk[5].strip()) == 0 and len(modelparams.chem_sch_mrk[2].strip()) == 0):
		# 		pr_flag = 1
					
		# if (pr_flag == 1):
		# 	# get the elements in line separated by peroxy radical separator
		# 	line2 = line1.split(modelparams.chem_sch_mrk[2])
		# 	RO2_count += 1 # count on number of lines considered in peroxy radical list
			
		# 	for line3 in line2: # loop through elements in line
		# 		if len(line3.split('='))>1: # in case of RO2 = ...
		# 			line3 = (line3.split('='))[1]
		# 		if len(line3.split(';'))>1: # in case of RO2 list finishing with ...;
		# 			line3 = (line3.split(';'))[0]
		# 		if len(line3.split('&'))>1: # in case of RO2 list finishing with &
		# 			line3 = (line3.split('&'))[0]
				
		# 		# remove any white space
		# 		line3 = line3.strip()
		# 		# don't include white space or ampersands
		# 		if (line3 == '' or line3 == '&'):
		# 			continue
				
		# 		else: # if not these, then first strip surrounding marks
		# 			if line3[0:len(modelparams.chem_sch_mrk[3])] == modelparams.chem_sch_mrk[3]:
		# 				line3 = line3[len(modelparams.chem_sch_mrk[3])::]
		# 			if line3[-len(modelparams.chem_sch_mrk[4])::] == modelparams.chem_sch_mrk[4]:
		# 				line3 = line3[0:-len(modelparams.chem_sch_mrk[4])]
				
		# 			RO2_names.append(line3)
		# 	# check for end of RO2 list - given either by marker for end or absence of
		# 	# marker for continuation onto next line of RO2
		# 	# check for marker for end of RO2 list
		# 	if (len(modelparams.chem_sch_mrk[5].strip()) > 0):
		# 		# .* allows search across all elements of line, \\ ensures marker is
		# 		# recognised as string
		# 		mark = str('.*\\' + modelparams.chem_sch_mrk[5])
		# 		if (re.match(mark, line1) != None):
		# 			pr_flag = 0
			
		# 	else: # look for absence of marker for RO2 list continuing onto next line
		# 		# .* allows search across all elements of line, \\ ensures marker is
		# 		# recognised as string
		# 		mark = str('.*\\' + modelparams.chem_sch_mrk[6])
		# 		if (re.match(mark, line1) == None):
		# 			pr_flag = 0
		
		# --------------------------------------------------------------------------------
		# gas-phase reaction equation part
		# ^ means occurs at start of line and, first \ means second \ can be interpreted 
		# and second \ ensures recognition of marker
		marker = str('^\\' +  modelparams.chem_sch_mrk[-1])
		
		# first check is whether equation start marker is present
		if (re.match(marker, line1) != None):
			
			# second check is whether markers for starting reaction rate coefficients
			# part, and markers for end of equation lines, are present
			eqn_markers = [str('.*\\' +  modelparams.chem_sch_mrk[9]), str('.*\\' +  modelparams.chem_sch_mrk[11])]
			if (re.match(eqn_markers[0], line1) != None and 
				re.match(eqn_markers[1], line1) != None):
				modelparams.eqn_list.append(line1) # store reaction equations

		# aqueous-phase reaction equation part
		# ^ means occurs at start of line and, first \ means second \ can be interpreted 
		# and second \ ensures recognition of marker
		# first, check if a marker given, if not bypass
		if modelparams.chem_sch_mrk[8] == '':
			continue
		else:
			
			marker = str('^\\' +  modelparams.chem_sch_mrk[8])
		
			if (re.match(marker, line1) != None):
				# second check is whether markers for starting reaction rate coefficients
				# part, and markers for end of equation lines, are present
				eqn_markers = [str('.*\\' +  modelparams.chem_sch_mrk[9]), str('.*\\' +  modelparams.chem_sch_mrk[11])]
				if (re.match(eqn_markers[0], line1) != None and 
					re.match(eqn_markers[1], line1) != None):
					modelparams.aqeqn_list.append(line1) # store reaction equations

		# surface (e.g. wall) reaction equation part
		# ^ means occurs at start of line and, first \ means second \ can be interpreted 
		# and second \ ensures recognition of marker
		# first, check if a marker given, if not bypass
		if (modelparams.chem_sch_mrk[12] == ''):
			continue
		else:
			
			marker = str('^\\' +  modelparams.chem_sch_mrk[12])
		
			if (re.match(marker, line1) != None):
				# second check is whether markers for starting reaction rate coefficients
				# part, and markers for end of equation lines, are present
				eqn_markers = [str('.*\\' +  modelparams.chem_sch_mrk[9]), str('.*\\' +  modelparams.chem_sch_mrk[11])]
				if (re.match(eqn_markers[0], line1) != None and 
					re.match(eqn_markers[1], line1) != None):
					modelparams.sueqn_list.append(line1) # store reaction equations

	# number of equations
	modelparams.eqn_num = np.array((len(modelparams.eqn_list), len(modelparams.aqeqn_list), len(modelparams.sueqn_list)))
    
	return(rrc, rrc_name, RO2_names, modelparams)
