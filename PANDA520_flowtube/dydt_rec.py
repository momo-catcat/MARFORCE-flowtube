'''module for calculating and recording change tendency (# molecules/cm3/s) of components'''
# changes due to gas-phase photochemistry and partitioning are included; 
# generated in init_conc and treats loss from gas-phase as negative

# File Created at 2022-08-03 10:54:09.539053

import numpy as np 

def dydt_rec(y, rindx, rstoi, reac_coef, pindx, pstoi, nprod, step, dydt_vst, nreac, num_sb, num_comp, pconc, core_diss, Psat, kelv_fac, kimt, kw, Cw, act_coeff, seedi, dydt_erh_flag, H2Oi, wat_hist):
	
	# loop through components to record the tendency of change, note that components can be grouped, e.g. RO2 for non-HOM-RO2 
	dydtnames = dydt_vst['comp_names'] 
	for comp_name in dydtnames: # get name of this component
		
		key_name = str(str(comp_name) + '_comp_indx') # get index of this component
		compi = dydt_vst[key_name] 
		# open relevant dictionary value containing reaction numbers and results 
		key_name = str(str(comp_name) + '_res')
		dydt_rec = dydt_vst[key_name] 
		# open relevant dictionary value containing flag for whether component is reactant (1) or product (0) in each reaction 
		key_name = str(str(comp_name) + '_reac_sign')
		reac_sign = dydt_vst[key_name] 
		# keep count on relevant reactions 
		reac_count = 0 
		# loop through relevant reactions 
		for i in dydt_rec[0, 0:-2]: # final two rows for particle- and wall-partitioning 
			i = int(i) # ensure reaction index is integer - this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float) 
			# estimate gas-phase change tendency for each reaction involving this component 
			gprate = ((y[rindx[i, 0:nreac[i]]]**rstoi[i, 0:nreac[i]]).prod())*reac_coef[i] 
			dydt_rec[step+1, reac_count] += reac_sign[reac_count]*((gprate))
			reac_count += 1 # keep count on relevant reactions 
			
		# now estimate and record tendency to change due to particle- and wall-partitioning  
		# particle-partitioning 
		# if efflorescence has occurred (modelled inside act_coeff_update), then need
		# to account for the immediate transfer of water from particles to gas
		if (wat_hist == 0):
			if (dydt_erh_flag != 0): 
				dydt_rec[step+1, reac_count] += dydt_erh_flag
			else: 
				dydt_rec[step+1, reac_count] = 0
		else:
			for ibin in range(num_sb-1): # size bin loop
				Csit = y[num_comp*(ibin+1):num_comp*(ibin+2)]
				conc_sum = np.zeros((1)) 
				if (sum(sum(pconc > 0.)) > 0): # if seed particles present 
					conc_sum[0] = ((Csit.sum()-sum(Csit[seedi]))+sum(Csit[seedi]*core_diss))
				else: 
					conc_sum[0] = Csit.sum() 
				# prevent numerical error due to division by zero 
				ish = (conc_sum == 0.) 
				conc_sum[ish] = 1.e-40 
				# particle surface gas-phase concentration (molecules/cc (air)) 
				Csit = (Csit[compi]/conc_sum)*Psat[0, compi]*kelv_fac[ibin, 0]*act_coeff[0, compi] 
				# partitioning rate (# molecules/cc.s) 
				dydt_all = kimt[ibin, compi]*(y[compi]-Csit) 
				# gas-phase change (# molecules/cm3/s) 
				dydt_rec[step+1, reac_count] -= np.sum(dydt_all) 
				
		# wall-partitioning 
		if (any(kw > 1.e-10)): 
			# concentration at wall (# molecules/cm3 (air)) 
			Csit = y[num_comp*num_sb:num_comp*(num_sb+1)] 
			Csit = (Psat[0, :]*(Csit/Cw)*act_coeff[0, :])
			dydt_all = (kw)*(y[0:num_comp]-Csit) 
			# gas-phase change (# molecules/cm3/s) 
			dydt_rec[step+1, reac_count+1] -= np.sum(dydt_all[compi]) 
		
	return(dydt_vst,dydt_rec) 
