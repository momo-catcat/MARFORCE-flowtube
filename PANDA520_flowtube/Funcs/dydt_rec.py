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
'''module for calculating and recording change tendency (# molecules/cm3/s)
 of components'''
# changes due to gas-phase photochemistry and partitioning are included; 
# generated in init_conc and treats loss from gas-phase as negative

# File Created at 2026-07-06 14:19:35.158753

import numpy as np 

def dydt_rec(y, reac_coef, step, modelparams):
	
	# loop through components to record the tendency of
	#change, note that components can be grouped, e.g. RO2 for non-HOM-RO2 
	dydtnames = modelparams.dydt_vst['comp_names'] 
	for comp_name in dydtnames: # get name of this component
		
		key_name = str(str(comp_name) + '_comp_indx') # get index of this component
		compi = modelparams.dydt_vst[key_name] 
		# open relevant dictionary value containing reaction numbers and results 
		key_name = str(str(comp_name) + '_res')
		dydt_rec = modelparams.dydt_vst[key_name]
		if hasattr(modelparams, 'sim_ci_file'):
			# open relevant dictionary value containing estimated continuous influx
			key_name = str(str(comp_name) + '_ci')
			ci_array = modelparams.dydt_vst[key_name]
		# open relevant dictionary value containing flag
		# for whether component is reactant (1) or product (0) in each reaction 
		key_name = str(str(comp_name) + '_reac_sign')
		reac_sign = modelparams.dydt_vst[key_name] 
		# keep count on relevant reactions 
		reac_count = 0 
		# loop through relevant reactions 
		# note that final three rows are for
		# particle- and wall-partitioning and dilution 
		for i in dydt_rec[0, 0:-3]: 
			i = int(i) # ensure reaction index is integer - this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float) 
			# estimate gas-phase change tendency for each reaction involving this component 
			gprate = ((y[modelparams.rindx_g[i, 0:modelparams.nreac_g[i]]]**modelparams.rstoi_g[i, 0:modelparams.nreac_g[i]]).prod())*reac_coef[i] 
			dydt_rec[step+1, reac_count] += reac_sign[reac_count]*((gprate))
			reac_count += 1 # keep count on relevant reactions 
			
		
		# dilution
		if (modelparams.dil_fac_now > 0):
			dydt_rec[step+1, reac_count+2] -= y[compi]*modelparams.dil_fac_now
	
	return(modelparams) 
