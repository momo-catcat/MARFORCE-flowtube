def cmd_calib5(const_comp_conc, params, Init_comp_conc, Q1, Q2):
    # %% import packages
    import numpy as np
    import RO2_conc
    import sch_interr
    import eqn_interr
    import eqn_pars
    import init_conc
    import RO2_indices
    import rate_coeffs
    import write_rate_file
    import cal_const_comp_conc
    from judg_spe_reac_rates import jude_species as jude_species
    from get_diff_and_u import get_diff_and_u_for_more_species
    from get_formula import get_formula
    from model_1 import model_1
    from model_3 import model_3
    from model_4 import model_4
    from model_5_test import model_5
    from meanconc_cal import meanconc_cal
    from grid_parameters import grid_para as grid_para

    # % load the inputs
    T = params['T']  # constant species
    p = params['p']  # temperature
    R1 = params['R1']  # diameters for first tube
    R2 = params['R2']  # diameters for second tube
    L1 = params['L1']  # length for first tube
    L2 = params['L2']  # length for second tube
    Q1 = Q1 / 60  # flow for first tube
    Q2 = Q2  / 60  # flow for second tube
    sch_name = params['sch_name']  # file for the MCM file
    chm_sch_mrk = ['{', 'RO2', '+', '', '', ';', '+', ';', '$', '{', ':', ';', '}'] #params['chm_sch_mrk']  # markers to isolate sections of chemical scheme based on MCM KPP format
    const_comp = params['const_comp']
    Init_comp = params['Init_comp']
    Diff_setname = params['Diff_setname']
    Diff_set = params['Diff_set']
    con_infl_nam = const_comp
    Zgrid = params['Zgrid']  # number of grid points in tube length direction
    Rgrid = params['Rgrid']  # number of grid points in tube radius direction
    key_spe_for_plot = params['key_spe_for_plot']  # key species for plot
    plot_spec = params['plot_spec']  # plot species
    dt = params['dt']
    flag_tube = params['flag_tube']
    const_comp_free = params['const_comp_free']  # the species in dilution flow
    const_comp_conc_free = params['const_comp_conc_free']

    # read the file and separate the equations and rate coefficients
    eqn_list, num_eqn, rrc, rrc_name, RO2_names, eqn_list_on = sch_interr.sch_interr(chm_sch_mrk, sch_name)
    #print(eqn_list_on)
    # find the comp_namelist reac_coef, and indx for products and reactants
    [rindx, rstoi, pindx, pstoi, reac_coef, nreac, nprod, comp_namelist, comp_num] = eqn_interr.eqn_interr(
        num_eqn, eqn_list, chm_sch_mrk)
    # find RO2 and the constant concentration
    RO2_indx, HOMRO2_indx, con_C_indx = eqn_pars.extr_mech(sch_name, chm_sch_mrk,
                                                           con_infl_nam, const_comp)

    RO2_indi = RO2_indices.RO2_indices(comp_namelist, RO2_names)
    # get the diffusion for all species and  the index of species in C except constant compounds
    print('comp',comp_namelist)
    print('Diff_set',Diff_set)
    print('Diff_setname', Diff_setname)
    u, Diff_vals = get_diff_and_u_for_more_species(comp_namelist, Diff_setname, con_C_indx, Diff_set, T, p)

    numLoop = 500  # number of times to run to reach the pinhole of the instrument
    timesteps = 1000 # number of timesteps, dt * timesteps * numLoop is time elapsed in the final solution

    # Change odd number Rgrid to even number grid
    if (Rgrid % 2) != 0:
        Rgrid = Rgrid + 1
    # % apply all the concentration for the constant comp
    const_comp_gird = cal_const_comp_conc.cal_const_comp_conc(Rgrid, Zgrid, const_comp_conc, L1, L2, const_comp)
    # % set the concentration for all the species for the grid of 80*40 in c
    c = np.zeros([Rgrid, Zgrid, comp_num])

    for i in Init_comp:  # set [OH] at z = 0 # set [HO2] at z = 0 [oh]
        c[:, 0, comp_namelist.index(i)] = Init_comp_conc[Init_comp.index(i)]
    for i in const_comp:  # set the constant concentrations for const_comp
        c[:, :, comp_namelist.index(i)] = const_comp_gird[const_comp.index(i)]

    # % write the rate coefficients in a new file (rate_coeffs.py)
    write_rate_file.write_rate_file(reac_coef, p, rrc, rrc_name, 0)
    # % store the first column of the initial concentration into C0
    C0 = c[0, 0, :]
    # % store indx for products and reactants to dydt_vst, and the initial concentration for first column
    [y, M, dydt_vst] = init_conc.init_conc(comp_num, comp_namelist, C0, T, p, comp_namelist, \
                                           rindx, pindx, num_eqn[0], nreac, nprod,
                                           comp_namelist, \
                                           RO2_indx, HOMRO2_indx, rstoi, pstoi)
    # % calculate the RO2 concentration
    RO2conc = RO2_conc.RO2_conc(RO2_indi, y)
    # % calculate H2O, O2, NO, HO2, NO3
    op = jude_species(y, comp_namelist)
    # % calculate reaction rate coefficients values

    rate_values = rate_coeffs.evaluate_rates(RO2conc, T, 0, M, M * 0.7809, op[0], op[1], op[2], op[3],
                                             op[4], p)
    # used as the title for the plotted figures
    formula = get_formula(plot_spec)
    # % set the grids parameters
    Rtot, dr, dx, sp_line = grid_para(Zgrid, Rgrid, R2, R1, L2, L1, comp_num)
    #print(rate_values)
    #print(comp_namelist)
    #print(comp_num)
    # %% run the modules and plot
    if flag_tube == '3':
        c = model_3(R2, R1, Rgrid, Zgrid, comp_num, L2, L1, numLoop, comp_namelist, key_spe_for_plot, dt, timesteps,
                    Diff_vals, Rtot, Q1, Q2, dydt_vst, rindx, nreac, rstoi, rate_values, const_comp, u, plot_spec,
                    formula, c, dr, dx,params['model_mode'],const_comp_conc)
    elif flag_tube == '4':
        c = model_4(R2, R1, Rgrid, Zgrid, comp_num, L2, L1, numLoop, comp_namelist, key_spe_for_plot, dt, timesteps,
                    Diff_vals, Rtot, const_comp_free, const_comp_conc_free, Q1, Q2, dydt_vst, rindx, nreac, rstoi,
                    rate_values,
                    const_comp, u, plot_spec, formula, c, dr, dx,params['model_mode'])
    else:  # model 1 and model 2 they are same
        # % run once with two different tubes or one tube
        c = model_1(R2, Rgrid, Zgrid, L2, L1, numLoop, comp_namelist, key_spe_for_plot, dt,
                    timesteps, Diff_vals, Rtot, Q1, dydt_vst, rindx, nreac, rstoi, rate_values, const_comp, u,
                    plot_spec, formula, c, dr, dx, params['model_mode'])
    # % calculate the meanconc for each species
    meanConc = meanconc_cal(R2, Rgrid, plot_spec, comp_namelist, c, params['model_mode'])

    return meanConc, c