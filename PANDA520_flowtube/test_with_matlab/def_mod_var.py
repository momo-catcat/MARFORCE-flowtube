'''default model variables'''
# the default model variables, loaded by PyCHAM on starting

import numpy as np
import os
import pickle # used for storing values

def def_mod_var(caller): # define function

	# inputs: -----------------------------------------------
	# caller - mark for calling function
	# -------------------------------------------------------

	# default input files ------------------------------
	# default chemical scheme
	# sch_name = 'C:/Users/jiali/PANDA520-flowtube/PANDA520_flowtube/input_mechanism/SO2_SA.txt'
	sch_name = '/Users/momo/Documents/science/coding/git/Projects/PANDA520-flowtube/PANDA520_flowtube/input_mechanism/SO2_SA.txt'
	# xml_name = 'C:/Users/jiali/PANDA520-flowtube/PANDA520_flowtube/input_mechanism/ex_xml_SO2.xml' # xml file path
	xml_name = '/Users/momo/Documents/science/coding/git/Projects/PANDA520-flowtube/PANDA520_flowtube/input_mechanism/ex_xml_SO2.xml'  # xml file path
	inname = 'Default' # model variables file name
	
	# general ---------------------------------------------------------------------------------
	# name of folder to save results to
	sav_nam = 'default_res_name'
	# markers to isolate sections of chemical scheme based on MCM KPP format
# 	chem_sch_mark = ['{', 'RO2', '+', 'C(ind_', ')','' , '&', '' , '', ':', '}', ';','}']
	chm_sch_mrk = ['{', 'RO2', '+', '', '', ';', '+', ';', '$', '{', ':', ';','}']
    
	# time interval between updates to integration inputs (s)
	update_stp = 1.
	tot_time = 1. # total time to integrate over (s)
	save_step = 1. # time interval between saving results (s)
	if (caller == 0): # called from PyCHAM
		uman_up = 0 # marker for whether to update the UManSysProp folder
	if (caller == 1): # called from Travis
		uman_up = 1
	int_tol = [1.e-3, 1.e-4] # integration tolerances (absolute first, relative second)
	testf = 0 # whether in testing mode or not
	
	# chamber environment -----------------------------------------------------------------
	temp = np.array((298.15)).reshape(1) # temperature of experiment (K)
	tempt = np.array((0.0)).reshape(1) # time that temperatures reached (s)	
	RH = np.array(([0.65])) # humidity of experiment (fraction of 1)
	RHt = np.array(([0])) # time through simulation (s) RH reached
	Press = 101300 # air pressure during experiment (Pa)
	dil_fac = 0. # dilution factor (volume fraction per second)
	drh_str = str('0.*TEMP')
	erh_str = str('0.*TEMP')
# 	# particle section ----------------------------------------------------------------
# 	siz_stru = 0 # size structure (0 for moving-centre, 1 for full-moving)
# 	num_sb = 0 # number of particle size bins
# 	# whether particle number concentrations expressed by modes (0) or explicitly	
# 	pmode = 0
# 	# concentration of particles (# particles/cc (air))
# 	pconc = np.zeros((1, 1))
# 	pconct = np.zeros((1, 1)) # times of particle injection (s)
# 	pcont = np.zeros((1, 1), dtype=int) # whether particle injection instantaneous or continuous
# 	seed_mw = (np.array((132.14))).reshape(1) # molecular weight of seed material (g/mol)
# 	seed_diss = [1.] # dissociation constant of seed material
# 	seed_dens = np.ones((1)) # density of seed material (g/cc)
# 	seed_name = ['core'] # name of component forming seed material
# 	seedx = np.ones((1, 1)) # mole fraction of dry seed components
# 	Vwat_inc = 1 # seed particle number size distribution does include water volume
# 	seed_eq_wat = 1 # allow water equilibration with seed particle before experiment starts
# 	lowsize = 0. # smallest size bin boundary (radius) (um)
# 	uppsize = 5.e-1 # largest size bin boundary (radius) (um)
# 	space_mode = 'lin' # treatment for spacing between size bins
# 	std = np.ones((1, 1))*1.2 # standard deviation for particle number size distribution
# 	mean_rad = np.ones((1, 1))*-1.e6 # mean radius for particle number size distribution (um)
# 	new_partr = 2.e-7 # radius of newly nucleated particles (cm)
# 	# nucleation parameters
# 	nucv1 = 0.
# 	nucv2 = 0.
# 	nucv3 = 0.
# 	nuc_comp = ['core'] # chemical scheme name of nucleating component
# 	# marker to say whether or not to adapt integration time interval 
# 	# and initial condition update to nucleation
# 	nuc_ad = 1
# 	ser_H2O = 1 # whether to serialise gas-particle partitioning of water
# 	coag_on = 1 # whether to model coagulation
# 	# flag for particle-phase history with respect to water (0 for dry and therefore 
# 	# on the deliquescence curve, 1 for wet and therefore on the efflorescence curve)
# 	wat_hist = 1
# 	# string describing the deliquescence relative humidity (fraction 0-1) dependence on 
# 	# temperature
 	# drh_str = str('0.*TEMP')
# 	# string describing the efflorescence relative humidity (fraction 0-1) dependence on 
# 	# temperature
# 	erh_str = str('0.*TEMP')
# 	# fraction of total gas-particle partitioning coefficient below which the 
# 	# partitioning coefficient is set to zero, e.g. because surface area of a size bin
# 	# is relatively very small
# 	z_prt_coeff = 1.e-9

	# gas inputs ------------------------------------------------------------
	# chemical scheme name of components present initially
	comp0 = np.array(())
	# initial concentrations (ppb)
	y0 = np.array(())	
	con_infl_nam = [] # chemical scheme names of components with continuous influx
	# influx rate of components with continuous influx (ppb/s)
	con_infl_C = np.empty(0)
	# times of component influx (s)
	con_infl_t = np.empty(0)
	# chemical scheme name of components with constant concentration	
	const_comp = []
	# Chemical scheme names of components injected instantaneously after start of experiment
	Compt = []
	# times at which instantaneous injection of component(s) occur after 
	# experiment start (s)
	injectt = np.empty(0)
	# concentration(s) (ppb) of component(s) injected instantaneously after 
	# experiment start
	Ct = np.zeros((0, 0))
	# the gas-particle partitioning cutoff (Pa)
	partit_cutoff = []

	# lights -------------------------------------------------------------------------------
# 	light_stat = np.zeros((1), dtype='int') # light status
# 	light_time = np.zeros((1, 1)) # time that light status attained (s)
# 	daytime = 0. # time of day experiment starts (s)
# 	lat = 0. # latitude of experiment (degrees)
# 	lon = 0. # longitude of experiment (degrees)
# 	af_path = 'no' # path to customised (non-MCM) actinic flux file
# 	# path to file containing absorption cross-sections and quantum yields
# # 	photo_path = str(os.getcwd() + '/PyCHAM/photofiles/' + 'MCMv3.2')
# 	dayOfYear = 1 # number of days since 31st December that experiment conducted
# 	tf = 1. # transmission factor for natural sunlight (0-1 fraction)
# 	# marker to say whether or not to adapt integration time interval 
# 	# and initial condition update to changing natural light intensity
# 	light_ad = 1
# 	tf_UVC = 1. # transmission factor for 254 nm wavelength artificial light

	# deposition of particles and vapours to wall ------------------------------------------
# 	wall_on = 1 # marker for whether to consider wall (0 for no, 1 for yes)
# 	Cw = 0. # effective absorbing mass of wall (g/m3 (air))
# 	kw = np.zeros((1)) # gas-wall mass transfer coefficient (/s)

# 	inflectDp = 0. # diameter of deposition function inflection
# 	pwl_xpre = 0. # gradient before inflection
# 	pwl_xpro = 0. # gradient after inflection
# 	inflectk = 0. # rate at inflection
# 	chamSA = 42. # chamber surface area (m2)
# 	chamV = 18. # chamber volume (m3, set to MAC value)
# 	Rader = -1 # flag for deposition to wall treatment (0 for customised, 1 for Rader and McMurry (1985))
# 	p_char = 0. # average number of charges per particle (/particle)
# 	e_field = 0. # average electric field inside chamber (g.m/A.s3)

	# specific component properties ---------------------------------------------------------
	# chemical scheme name of components to track the change tendencies of	
	dydt_trak = []
	# chemical scheme names of components with densities manually set
	dens_comp = []
	# manually assigned densities (g/cc)
	dens = []
	# chemical scheme names of components with vapour pressures manually set
	vol_comp = []
	# manually assigned vapour pressures (Pa)
	volP = []
	# names of components (corresponding to chemical scheme name) with 
	# activity coefficient stated in act_user
	act_comp = []
	# user-specified activity coefficients of components with names given 
	# in act_comp
	act_user = []
	# names of components with user-defined accommodation coefficients	
	accom_comp = []
	# user-defined accommodation coefficients
	accom_val = []
	

	# prepare for pickling
	list_vars = [sav_nam, sch_name, chm_sch_mrk, xml_name, inname, update_stp, 
			tot_time, comp0, y0, temp, tempt, RH, RHt, Press, 
            save_step, const_comp, Compt, injectt, Ct,con_infl_nam, con_infl_t, 
            con_infl_C, dydt_trak,
            dens_comp, dens, vol_comp, 
			volP, act_comp, act_user, accom_comp, accom_val, uman_up, 
			int_tol,dil_fac, partit_cutoff,drh_str, erh_str, testf]
# 	list_vars = [sav_nam, sch_name, chem_sch_mark, xml_name, inname, update_stp, 
# 			tot_time, comp0, y0, temp, tempt, RH, RHt, Press, wall_on, 
# 			Cw, kw, siz_stru, num_sb, pmode, pconc, pconct, lowsize, 
# 			uppsize, space_mode, std, mean_rad, save_step, const_comp, 
# 			Compt, injectt, Ct, seed_name, seed_mw, seed_diss, seed_dens, 
# 			seedx, light_stat, light_time, daytime, lat, lon, af_path, 
# 			dayOfYear, tf, light_ad, con_infl_nam, 
# 			con_infl_t, con_infl_C, dydt_trak, dens_comp, dens, vol_comp, 
# 			volP, act_comp, act_user, accom_comp, accom_val, uman_up, 
# 			int_tol, new_partr, nucv1, nucv2, nucv3, nuc_comp, nuc_ad, 
# 			coag_on, inflectDp, pwl_xpre, pwl_xpro, inflectk, chamSA, 
# 			Rader, p_char, e_field, dil_fac, partit_cutoff, ser_H2O, 
# 			wat_hist, drh_str, erh_str, pcont, Vwat_inc, seed_eq_wat, 
# 			z_prt_coeff, tf_UVC, testf, chamV]

	
	# path to store for variables
	input_by_sim = str(os.getcwd() + '/pickle.pkl')
		
	with open(input_by_sim, 'wb') as f: # the file to be used for pickling
		pickle.dump(list_vars,f) # pickle
		f.close() # close


	return(sav_nam, sch_name, chm_sch_mrk, xml_name, inname, update_stp, 
		tot_time, comp0, y0, temp, tempt, RH, RHt, Press,  
		save_step, const_comp, Compt, injectt, Ct,  
		con_infl_nam, con_infl_t, con_infl_C, 
		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
		accom_comp, accom_val, uman_up, int_tol, dil_fac, partit_cutoff,drh_str, erh_str, testf)
