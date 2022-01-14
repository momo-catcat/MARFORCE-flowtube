def cmd_calib5(const_comp_conc, R1, L1, Q1, R2, L2, Q2, It, fullOrSimpleModel, time):       
    #%% import packages
    
    import numpy as np
    import os 
    import pickle
    import RO2_conc
    import rate_coeffs
    import sch_interr
    import eqn_interr
    import eqn_pars
    import water_calc
    import init_conc
    import RO2_indices  
    import write_rate_file
    import cal_const_comp_conc
    from Vapour_calc import H2O_conc
    import matplotlib.pyplot as plt
    from scipy import interpolate
    from odesolve3 import odesolve as odesolve

    #%% sort chemistry   
    
    # load the pick file for the inital inputs 
    input_by_sim = os.getcwd()+'/pickle.pkl'

    with open(input_by_sim, 'rb') as pk:
        
        [sav_nam, sch_name, chm_sch_mrk, xml_name, 
    			comp0, y0, T,  p, const_comp,  con_infl_nam,  
                dydt_trak, uman_up, drh_str, erh_str, testf]= pickle.load(pk) 
        
    pk.close()
    
	# read the file and store everything into a list   
    f_open_eqn = open(sch_name, mode='r') # open the chemical scheme file

    total_list_eqn = f_open_eqn.readlines()
    f_open_eqn.close() # close file
    
   
    eqn_list, aqeqn_list, num_eqn, rrc, rrc_name, RO2_names, eqn_list_on=sch_interr.sch_interr(total_list_eqn, chm_sch_mrk)
    
    [rindx, rstoi, pindx, pstoi, reac_coef, 
			nreac, nprod, y_arr, y_rind, uni_y_rind, y_pind, 
			uni_y_pind, reac_col, prod_col, rstoi_flat, pstoi_flat, 
			rr_arr, rr_arr_p, comp_namelist, comp_list, Pybel_objects, 
			comp_num] = eqn_interr.eqn_interr(num_eqn, eqn_list, aqeqn_list, chm_sch_mrk)
    
    # comp_namelist, comp_list
    # nreac :  number of the reactant 
    # nprod : number of production A + B = C nprod is 1 nreac is 2 
    # rstoi_flat : record 1D array of stoichiometries per equation A + B = C is 1 1 1,     
    [rindx, pindx, rstoi, pstoi, nreac, nprod, y_arr, y_rind,
    		uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col, 
    		rstoi_flat, pstoi_flat, rr_arr, rr_arr_p, 
    		comp_num, RO2_indx, 
    		HOMRO2_indx, comp_list, 
    		Pybel_objects, eqn_num, comp_namelist, 
    		erf, err_mess, con_C_indx] = eqn_pars.extr_mech(sch_name, chm_sch_mrk, xml_name, 
       const_comp, drh_str, erh_str,  sav_nam)
                                                                    
    dydt_trak=comp_namelist  
                                                                                  
    RO2_indices = RO2_indices.RO2_indices(comp_namelist, RO2_names)

    const_comp_index = [comp_namelist.index(const_comp[i]) for i in range(len(const_comp))]
    
    u = list(range(len(comp_namelist))) #% u is the index of species in C except constant compounds 
    for i in reversed(const_comp_index):
        del u[i]                                                                     
              

    H2Oconc = const_comp_conc[:,const_comp.index('H2O')]                                                       
    # if H2Oconc[0] == 0:
    #     meanWeightedH2SO4 = 0
        
    if H2Oconc[1] == 0:
        H2Oconc[1] = H2Oconc[0]*Q1/Q2
    
    # Initial OH concentration
    csH2O = 7.22e-20 #cm2
    
    qyH2O = 1
    #It=1.84e10
    OHconc = It * csH2O * qyH2O * H2Oconc[0]

    # diffusion constants, [cm**2/s]
    rh =  H2Oconc[0] * 1e6 * 1.3806488e-23 * T / H2O_conc(T - 273.15, 1).SatP[0]
    # ['OH', 'SO2', 'HSO3', 'HO2', 'H2O', 'O2', 'H2O2', 'SO3', 'SA']
    D = np.array([0.215,0.1, 0.126, 0.141, 0.1, 0.126, 0.1,0.126, 0.08])
    # D = np.array([0.126, 0.126, 0.141, 0.08, 0.215]) #todo need to calculate Diffusion coefficient properly
    T0 = np.array([300, 300, 298, 298, 298,298,298,298,298])
    D = 101325 / p * D * ((T ** (3 / 2)) / (T0 ** (3 / 2)))


    H2O, Psat_water, H2O_mw = water_calc.water_calc(T, rh)

    dt = 0.0001                        # timestep [s]
    numLoop = 500                      # number of times to run to reach the pinhole of the instrument
    timesteps = 1000                    # number of timesteps, dt * timesteps * numLoop is time elapsed in the final solution

    Zgrid = np.array(40).astype(int)                         # number of grid points in tube length direction
    Rgrid = np.array(80)                       # number of grid points in tube radius direction

    #Change odd number Rgrid to even number grid
    if (Rgrid % 2) != 0:
        Rgrid = Rgrid + 1


    #% set the dr dx Q parameters for the tube 
    dr = np.zeros([int(Rgrid),int(Zgrid),comp_num])
    dx = (L1+L2) / (Zgrid - 1)
    dr[:,0:int(Zgrid*L1/(L2+L1)),:] = 2 * R1 / (Rgrid - 1)
    dr[:,int(Zgrid*L1/(L2+L1)):,:] =  2 * R2 / (Rgrid - 1)
    
    Qtot = np.zeros([int(Rgrid),int(Zgrid),comp_num])
    Qtot[:,0:int(Zgrid*L1/(L2+L1)),:] = Q1
    Qtot[:,int(Zgrid*L1/(L2+L1)):,:] =  Q2
    
    Rtot = np.zeros([int(Rgrid),int(Zgrid),comp_num])
    Rtot[:,0:int(Zgrid*L1/(L2+L1)),:] = R1
    Rtot[:,int(Zgrid*L1/(L2+L1)):,:] =  R2
    
    # make the const comp to one 
    const_comp_gird = cal_const_comp_conc.cal_const_comp_conc(Rgrid, Zgrid,const_comp_conc ,L1,L2,const_comp)
          
    #% initial conditions and set the const_comp conc

    c = np.zeros([Rgrid, Zgrid, comp_num])
    c[:, 0, 0] = OHconc  # set [OH] at z = 0
    c[:, 0, 3] = OHconc  # set [HO2] at z = 0. This equals OH conc
    for i in const_comp:
        c[:,:,comp_namelist.index(i)] = const_comp_gird[const_comp.index(i)]  


    write_rate_file.write_rate_file(reac_coef,p, rrc, rrc_name, 0)
    
    Comp0 = comp_namelist

    C0 = c[0,0,:] 

    [y,  y_mw, num_comp, M, y_indx_plot, dydt_vst, 
    comp_namelist,  erf, err_mess, NOi, HO2i, NO3i]=init_conc.init_conc(comp_num, Comp0, C0, T, C0[comp_namelist.index('H2O')], p, Pybel_objects,
    testf, dydt_trak, rindx, pindx, num_eqn[0], nreac, nprod, comp_namelist,  comp_namelist,  comp_namelist, RO2_indx, HOMRO2_indx, rstoi, pstoi)                                                                                              
    y=c[0,0,:]

    RO2conc = RO2_conc.RO2_conc(RO2_indices,y)
    rate_values, erf, err_mess = rate_coeffs.evaluate_rates(RO2conc, y[comp_namelist.index('H2O')], T, 0, M, M*0.7809, y[comp_namelist.index('O2')], 0,  y[comp_namelist.index('HO2')],0,p)
    formula = ['OH','$\mathdefault{HSO_3}$','$\mathdefault{HO_2}$','$\mathdefault{SO_3}$','SA']
    #%% plot

    for j in range(numLoop):
        c1 = c.copy()
        oldH2SO4 = c1[:, -1 ,comp_namelist.index('SA')]

        c = odesolve(timesteps, Zgrid, Rgrid, dt,  D, Rtot, dr, dx, Qtot,c,comp_namelist,dydt_vst,rindx,nreac,rstoi,rate_values,const_comp,u)

        t = ['OH','HSO3','HO2','SO3','SA']

        tim = (j + 1) * timesteps * dt
        comp_plot_index = [comp_namelist.index(t[i]) for i in range(len(t))]

        newH2SO4 = c[:, -1, comp_namelist.index('SA')]
        
        fig, axs = plt.subplots(2,3, figsize=(8, 5), facecolor='w', edgecolor='k')
        fig.subplots_adjust(hspace = .5, wspace=.45)
        plt.style.use('default')
        plt.rcParams.update({'font.size':13,'font.weight':'bold','font.family':'serif','font.serif':'Times New Roman'})

        axs = axs.ravel()
        
        for i in range(len(t)):
            axs[i].pcolor(np.linspace(0, L2+L1, Zgrid),np.linspace(-R1, R1, Rgrid), c[:, : , comp_plot_index[i]], shading = 'nearest')
            axs[i].pcolor(np.linspace(0, L2+L1, Zgrid),np.linspace(-R1, R2, Rgrid), c[:, : , comp_plot_index[i]], shading = 'nearest')

            axs[i].set_xlabel('L [cm]')
            axs[i].set_ylabel('R [cm]')
            axs[i].set_title(formula[i])

        fig.delaxes(axs[5])
        plt.gcf().text(0.7,0.3, 'Time = ' + str(tim) ,fontsize = 15)
        plt.draw()
        plt.pause(1)  

        
        # # print(newH2SO4)

        print(['t = ' + str(tim) + "  H2SO4 difference: " + str(np.sum(newH2SO4 - oldH2SO4))])
            
        
        if (j > 15) & (np.sum(newH2SO4 - oldH2SO4) / np.sum(oldH2SO4) < 1e-5):
            break

    dr_final = R2 / (Rgrid - 1) * 2
    x = np.arange(0, R2, dr_final) + dr_final
    y = np.flip(c[0 : int(Rgrid / 2), -1,comp_namelist.index('SA')])
    y1 = np.flip(c[0: int(Rgrid /2 ),-1, comp_namelist.index('HO2')])
    splineres = interpolate.splrep(x, y)
    splineres_HO2 = interpolate.splrep(x, y1)
    
    plt.subplot(2,3,6)
    plt.plot(np.arange(0, R2, dr_final), interpolate.splev(np.arange(0, R2, dr_final), splineres))
    plt.title('[$\mathdefault{H_2SO_4}$] at end of tube')
    plt.xlabel('r [cm]')
    plt.ylabel('Concentration, [cm$^{-3}$]')

    rVec = np.arange(0, R2, 0.001)
    cVec = interpolate.splev(rVec, splineres)
    cVec_HO2 = interpolate.splev(rVec, splineres_HO2)
    meanH2SO4 = 2 * 0.001 / R2 ** 2 * np.sum(cVec * rVec) #calculate average concentration over the cross section
    # meanWeightedH2SO4 = 4 * 0.001 / R2 ** 2 * np.sum(cVec * rVec * (1 - rVec ** 2 / R2 ** 2)) #not sure about the formulation. Needs to be checked
    
    meanHO2 = 2 * 0.001 / R2 ** 2 * np.sum(cVec_HO2 * rVec) #calculate average concentration over the cross section
    # meanWeightedHO2 = 4 * 0.001 / R2 ** 2 * np.sum(cVec_HO2 * rVec * (1 - rVec ** 2 / R2 ** 2)) #not sure about the formulation. Needs to be checked
    
    return(meanH2SO4,meanHO2,c)
