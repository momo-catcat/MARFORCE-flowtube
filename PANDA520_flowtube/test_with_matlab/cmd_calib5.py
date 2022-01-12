def cmd_calib5(O2conc, H2Oconc, SO2conc, R, L, Q,  It, T, p, fullOrSimpleModel, time):
    #%% parameters (that can be changed by the user)
    #import packages
    import numpy as np
    import os 
    import pickle
    import def_mod_var
    import chem_sch_SMILES
    import sch_interr
    import eqn_interr
    import eqn_pars
    # import scipy.constants as si
    # import rrc_calc  
    # import rate_coeffs 
    import water_calc
    import init_conc
    import RO2_indices  
    import write_rate_file
    from Vapour_calc import H2O_conc
    import matplotlib.pyplot as plt
    from scipy import interpolate
    from odesolve3 import odesolve as odesolve
    # R= 0.78
    # R1= 0.78
    # L= 55
    # L1= 68
    # Q=11*1000/60
    # Q1=11*1000/60
    R = 24/10/2 # mm the inner diameters of the tube
    L = 93 # mm
    Q = 11*1000/60 # lpm

    # save all the parameters in to the pickle file 
    def_mod_var.def_mod_var(0)

    # load the pick file 
    input_by_sim = os.getcwd()+'/pickle.pkl'

    with open(input_by_sim, 'rb') as pk:
       [sav_nam, sch_name, chm_sch_mrk, xml_name, inname, update_stp, 
   		tot_time, comp0, y0, temp, tempt, RH, RHt, Press,  
   		save_step, const_comp, Compt, injectt, Ct,  
   		con_infl_nam, con_infl_t, con_infl_C, 
   		dydt_trak, dens_comp, dens, vol_comp, volP, act_comp, act_user, 
   		accom_comp, accom_val, uman_up, int_tol, dil_fac, partit_cutoff,drh_str, erh_str, testf] = pickle.load(pk) 
    pk.close()
    
    f_open_eqn = open(sch_name, mode='r') # open the chemical scheme file
	# read the file and store everything into a list
    total_list_eqn = f_open_eqn.readlines()
    f_open_eqn.close() # close file
    const_comp = ['SO2','O2','H2O']
# 	comp_name = [] # list for chemical scheme names of components in the chemical scheme
# 	comp_smil = [] # list for the SMILE strings of components present in the chemical scheme
    comp_name, comp_smil, err_mess, H2Oi=chem_sch_SMILES.chem_scheme_SMILES_extr(sch_name, xml_name, chm_sch_mrk)
    
    eqn_list, aqeqn_list, num_eqn, rrc, rrc_name, RO2_names, eqn_list_on=sch_interr.sch_interr(total_list_eqn, chm_sch_mrk)
    
    [rindx, rstoi, pindx, pstoi, reac_coef, nreac, nprod, jac_stoi, jac_den_indx, njac_g, jac_indx, 				
    y_arr, y_rind, uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, rstoi_flat, pstoi_flat, 
    rr_arr, rr_arr_p,comp_namelist, comp_list, Pybel_objects,comp_num, RO_indx] = eqn_interr.eqn_interr(num_eqn, eqn_list, aqeqn_list, chm_sch_mrk, comp_name, comp_smil)

# comp_namelist, comp_list
# nreac :  number of the reactant 
# nprod : number of production A + B = C nprod is 1 nreac is 2 
# rstoi_flat : record 1D array of stoichiometries per equation A + B = C is 1 1 1, 

    [rindx, pindx, rstoi, pstoi, nreac, nprod, jac_stoi, 
		njac, jac_den_indx, jac_indx, y_arr, y_rind,
		uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, 
		rstoi_flat, pstoi_flat, rr_arr, rr_arr_p, 
		comp_num, RO2_indx, RO_indx,
		HOMRO2_indx, comp_list, 
		Pybel_objects, eqn_num, comp_namelist, 
		comp_name, comp_smil, erf, err_mess, con_C_indx] = eqn_pars.extr_mech(sch_name, chm_sch_mrk, xml_name, 
        comp_name, int_tol,  const_comp,
		drh_str, erh_str, dil_fac, sav_nam)
                                                                        
    dydt_trak=comp_namelist

                                                                                      
    RO2_indices = RO2_indices.RO2_indices(comp_namelist, RO2_names)                                                                                     

                                                                    
                                                                        
    if H2Oconc == 0:
        meanWeightedH2SO4 = 0 ## change the value
        
    # if H2Oconc_1 == 0:
    #     H2Oconc_1 = H2Oconc*Q/Q1
    # reaction constants (using cm**3 and s)
    kSO2pOH = 1.32e-12 * (T / 300) ** -0.7
    kOHpHO2 = 4.8e-11 * np.exp(250 / T)
    kOHpOH = 6.9e-31 * (T / 300) ** -0.8 * p / 1.3806488e-23 / T / 1e6
    kSO3p2H2O = 3.9e-41 * np.exp(6830.6 / T)
    kHSO3pO2 = 1.3e-12 * np.exp(-330 / T)
    
    # Initial OH concentration
    csH2O = 7.22e-20 #cm2
    qyH2O = 1
    #It=1.84e10
    OHconc = It * csH2O * qyH2O * H2Oconc

    # diffusion constants, [cm**2/s]
    # order is: HSO3, SO3, HO2, H2SO4, OH

    rh = H2Oconc * 1e6 * 1.3806488e-23 * T / H2O_conc(T - 273.15, 1).SatP[0]
    # ['OH', 'SO2', 'HSO3', 'HO2', 'H2O', 'O2', 'H2O2', 'SO3', 'SA']
    D = np.array([0.215,0.1, 0.126, 0.141, 0.1, 0.126, 0.1,0.126, 0.08])
    # ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']
    # D = np.array([0.126, 0.126, 0.141, 0.08, 0.215]) #todo need to calculate Diffusion coefficient properly
    # D = np.array([0.4, 0.4, 0.4, 0.4, 0.4])
    T0 = np.array([300, 300, 298, 298, 298,298,298,298,298])
    D = 101325 / p * D * ((T ** (3 / 2)) / (T0 ** (3 / 2)))
    # NA = si.Avogadro
    H2O, Psat_water, H2O_mw = water_calc.water_calc(T, rh)

    dt = 0.0001                        # timestep [s]
    numLoop = 500                      # number of times to run to reach the pinhole of the instrument
    timesteps = 1000                    # number of timesteps, dt * timesteps * numLoop is time elapsed in the final solution

    Zgrid = np.array(40).astype(int)                         # number of grid points in tube length direction
    Rgrid = np.array(80)                       # number of grid points in tube radius direction

    #Change odd number Rgrid to even number grid
    if (Rgrid % 2) != 0:
        Rgrid = Rgrid + 1

    #% # do not change
    # while np.gcd(Rgrid - 1, 10) != 1:  # make sure we don't get a grid point for r = 0 (would give Inf in calculation)
    #     Rgrid = int(Rgrid + 1)
    # comp_num  = 5
    # dr = np.zeros([int(Rgrid),int(Zgrid),comp_num])
    # dx = (L) / (Zgrid - 1)
    # dr[:,0:int(Zgrid*L/(L+L1)),:] = 2 * R / (Rgrid - 1)
    # dr[:,int(Zgrid*L/(L+L1)):,:] =  2 * R1 / (Rgrid - 1)
    
    # H2Otot = np.zeros([int(Rgrid),int(Zgrid)])
    # H2Otot[:,0:int(Zgrid*L/(L+L1))] = H2Oconc
    # H2Otot[:,int(Zgrid*L/(L+L1)):] =  H2Oconc_1
    
    # Qtot = np.zeros([int(Rgrid),int(Zgrid),comp_num])
    # Qtot[:,0:int(Zgrid*L/(L+L1)),:] = Q
    # Qtot[:,int(Zgrid*L/(L+L1)):,:] =  Q1
    
    # SO2tot = np.zeros([int(Rgrid),int(Zgrid)])
    # SO2tot[:,0:int(Rgrid*L/(L+L1))] = SO2conc
    # # SO2tot[:,int(Zgrid*L/(L+L1)):] =  SO2conc/2
    # SO2tot[:,int(Zgrid*L/(L+L1)):] =  SO2conc
    # O2tot = np.zeros([int(Rgrid),int(Zgrid)])
    # O2tot[:,0:int(Rgrid*L/(L+L1))] = O2conc
    # # O2tot[:,int(Zgrid*L/(L+L1)):] =  O2conc/2
    # O2tot[:,int(Zgrid*L/(L+L1)):] =  O2conc
    # Rtot = np.zeros([int(Rgrid),int(Zgrid),comp_num])
    # Rtot[:,0:int(Rgrid*L/(L+L1)),:] = R
    # Rtot[:,int(Zgrid*L/(L+L1)):,:] =  R1
    #% initial conditions
    # order in c is: HSO3, SO3, HO2, H2SO4, OH
    c = np.zeros([Rgrid, Zgrid, comp_num])
    c[:, 0, 0] = OHconc  # set [OH] at z = 0
    c[:, 0, 3] = OHconc  # set [HO2] at z = 0. This equals OH conc
    c[:,:,1] = SO2conc
    c[:,0,8] = OHconc
    c[:,:,5] = O2conc
    c[:,:,4] = H2Oconc
    
    # c[:,0:int(Rgrid*L/(L+L1)),:] = c [:,0:int(Rgrid*L/(L+L1)),:]
    ## check and output parameters to file
    # D_temp = np.tile(D, (Rgrid * Zgrid))
    # D = np.transpose(np.reshape(D_temp, (Rgrid, Zgrid, comp_num)), (0, 1, 2))

    # r = np.zeros([int(Rgrid), int(Zgrid),  comp_num])
    # r = np.abs(np.array([r[i, :, :] + i - r.shape[0] / 2 + 0.5 for i in range(r.shape[0])])) * dr
    write_rate_file.write_rate_file(reac_coef,p, rrc, rrc_name, 0)
    Comp0 = comp_namelist
    # C0 = [1,1,1,1,1,1,1,1,1]
    C0 = c[0,0,:] 
    # comp_num  = 9
    [y,  y_mw, num_comp, M, y_indx_plot, dydt_vst, 
    comp_namelist,  erf, err_mess, NOi, HO2i, NO3i]=init_conc.init_conc(comp_num, Comp0, C0, temp, C0[comp_namelist.index('H2O')], Press, Pybel_objects,
    testf, dydt_trak, rindx, pindx, num_eqn[0], nreac, nprod, comp_namelist, Compt, comp_namelist, comp_smil, comp_namelist, RO2_indx, HOMRO2_indx, rstoi, pstoi)                                                                                              
    y=c[0,0,:]
    import RO2_conc
    import rate_coeffs
    RO2conc = RO2_conc.RO2_conc(RO2_indices,y)
    rate_values, erf, err_mess = rate_coeffs.evaluate_rates(RO2conc, y[comp_namelist.index('H2O')], T, 0, M, M*0.7809, y[comp_namelist.index('O2')], 0,  y[comp_namelist.index('HO2')],0,p)

    #%% plt.figure(1, figsize = [8,6])
    for j in range(numLoop):

        c1 = c.copy()
        oldH2SO4 = c1[:, -1 ,comp_namelist.index('SA')]
        # print('oldH2SO4',oldH2SO4)


        c = odesolve(timesteps, Zgrid, Rgrid, dt, D, R, L, Q,  comp_namelist,dydt_vst,rindx,nreac,rstoi,rate_values,c)
       
        newH2SO4 = c[:, -1, comp_namelist.index('SA')]
        
        # print('newH2SO4',newH2SO4)
        # print('oldH2SO4',oldH2SO4)

        
        # t = ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']
        # ['OH', 'SO2', 'HSO3', 'HO2', 'H2O', 'O2', 'H2O2', 'SO3', 'SA']
        t = ['OH','HSO3','HO2','SO3','SA']
        # t= comp_namelist
        # for i in range(5):
        #     c[int(c[:, 0, i].size / 2) : -1, :, i] = np.flipud(c[1 : int(c[:, 0, i].size / 2), :, i]) # symmetry, program only returns values for 0..R

        tim = (j + 1) * timesteps * dt
        
        # y = np.linspace(-R, R, Rgrid)
        
        # y1 =Rtot[:,:,1]
        # y2 = y[1,:]
        # if tim in np.linspace(0, 48, 17):
            # ct = c[:, : ,i]
        order= [0,2,3,7,8]
        for i in range(5):
            # plt.figure(figsize=(10,8))
            plt.subplot(2,3, i + 1)
            plt.pcolor(np.linspace(0, L, Zgrid),np.linspace(-R, R, Rgrid), c[:, : ,order[i]], shading = 'nearest')
            # plt.pcolor(np.linspace(0, L+L1, Zgrid),np.linspace(-R1, R1, Rgrid), c[:, : , order[i]], shading = 'nearest')
            # plt.pcolor(c[:, : , i])
            #colorbar
            plt.xlabel('L [cm]')
            plt.ylabel('r [cm]')
            if i == 2:
                plt.title(['t =' + str(tim) + t[i]]) # also print time in subplot 2 (top-middle)
            else:
                plt.title(t[i])
        plt.draw()
        plt.pause(1)
        # dr_final = R1 / (Rgrid - 1) * 2
        # x = np.arange(0, R1, dr_final)
        # y = c[0 : int(Rgrid / 2), -1, 3]
        # y1 = c[0: int(Rgrid /2 ),-1, 2]
        # splineres = interpolate.splrep(x, y)
        # splineres_HO2 = interpolate.splrep(x, y1)
        # newH2SO4 = c[:, -1, 3]
        # plt.subplot(2,3,6)
        # plt.plot(np.arange(0, R1, dr_final), interpolate.splev(np.arange(0, R1, dr_final), splineres))
        # plt.title('[H2SO4] at end of tube')
        # plt.xlabel('r [cm]')
        # plt.ylabel('Concentration, [cm**{-3}]')
        newH2SO4 = c[:, -1, comp_namelist.index('SA')]

        # print('newH2SO4',newH2SO4)
        # print('oldH2SO4',oldH2SO4)
        # print('H2SO4',newH2SO4 - oldH2SO4)

        print(['time: ' + str(tim)])
        print(['t = ' + str(tim) + "  H2SO4 difference: " + str(np.sum(newH2SO4 - oldH2SO4))])
            
        
        if (j > 15) & (np.sum(newH2SO4 - oldH2SO4) / np.sum(oldH2SO4) < 1e-5):
            break

        dr_final = R / (Rgrid - 1) * 2
        x = np.arange(0, R, dr_final)
        y = c[0 : int(Rgrid / 2), -1,comp_namelist.index('SA')]
        y1 = c[0: int(Rgrid /2 ),-1, comp_namelist.index('HO2')]
        splineres = interpolate.splrep(x, y)
        splineres_HO2 = interpolate.splrep(x, y1)
        
        plt.subplot(2,3,6)
        plt.plot(np.arange(0, R, dr_final), interpolate.splev(np.arange(0, R, dr_final), splineres))
        plt.title('[H2SO4] at end of tube')
        plt.xlabel('r [cm]')
        plt.ylabel('Concentration, [cm**{-3}]')

        rVec = np.arange(0, R, 0.001)
        cVec = interpolate.splev(rVec, splineres)
        cVec_HO2 = interpolate.splev(rVec, splineres_HO2)
        meanH2SO4 = 2 * 0.001 / R ** 2 * np.sum(cVec * rVec) #calculate average concentration over the cross section
        meanWeightedH2SO4 = 4 * 0.001 / R ** 2 * np.sum(cVec * rVec * (1 - rVec ** 2 / R ** 2)) #not sure about the formulation. Needs to be checked
        
        meanHO2 = 2 * 0.001 / R ** 2 * np.sum(cVec_HO2 * rVec) #calculate average concentration over the cross section
        meanWeightedHO2 = 4 * 0.001 / R ** 2 * np.sum(cVec_HO2 * rVec * (1 - rVec ** 2 / R ** 2)) #not sure about the formulation. Needs to be checked
        
        # print(np.sum(cVec * rVec))
        # print(np. sum(cVec * rVec * (1 - rVec ** 2 / R1 ** 2)))

    return(meanH2SO4,meanHO2)
