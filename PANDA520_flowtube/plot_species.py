def plot(c,L1,L2,R1,R2,species):       
    import numpy as np
    import os 
    import pickle
    import sch_interr
    import eqn_interr
    import matplotlib.pyplot as plt
    # sort chemistry   
    
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
    
    Zgrid = np.array(40).astype(int)                         # number of grid points in tube length direction
    Rgrid = np.array(80).astype(int)                       # number of grid points in tube radius direction
    
    fig, axs = plt.subplots(1,1, figsize=(5,4), facecolor='w', edgecolor='k')
    plt.style.use('default')
    plt.rcParams.update({'font.size':13,'font.weight':'bold','font.family':'serif','font.serif':'Times New Roman'})
    axs.pcolor(np.linspace(0, L2+L1, Zgrid),np.linspace(-R1, R1, Rgrid), c[:, : , comp_namelist.index(species)], shading = 'nearest')
    axs.pcolor(np.linspace(0, L2+L1, Zgrid),np.linspace(-R1, R2, Rgrid), c[:, : , comp_namelist.index(species)], shading = 'nearest')
    # axs = axs.ravel()
    axs.set_xlabel('L [cm]')
    axs.set_ylabel('R [cm]')
    axs.set_title(species)

    plt.draw()
    plt.pause(1)  
    print (species+' last columns',c[:,-1,comp_namelist.index(species)])
   
    return(species)
