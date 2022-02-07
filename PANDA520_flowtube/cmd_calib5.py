def cmd_calib5(const_comp_conc, params, Init_comp_conc):
    # %% import packages
    import numpy as np
    import RO2_conc
    import sch_interr
    import eqn_interr
    import eqn_pars
    import init_conc
    import RO2_indices
    import write_rate_file
    import cal_const_comp_conc
    from judg_spe_reac_rates import jude_species as jude_species
    import matplotlib.pyplot as plt
    from scipy import interpolate
    from odesolve3_single import odesolve as odesolve_single
    from odesolve3 import odesolve as odesolve_double
    from odesolve3_Y import odesolve as odesolve_Y
    from get_diff_and_u import get_diff_and_u
    from get_formula import get_formula

    # % get the inputs
    T = params['T']  # constant species
    p = params['p']  # temperaure
    R1 = params['R1']  # diameters for frist tube
    R2 = params['R2']  # diameters for second tube
    L1 = params['L1']  # length for frist tube
    L2 = params['L2']  # length for second tube
    Q1 = params['Q1'] * 1000 / 60  # flow for frist tube
    Q2 = params['Q2'] * 1000 / 60  # flow for second tube
    sch_name = params['sch_name']  # file for the MCM file
    chm_sch_mrk = params['chm_sch_mrk']  # markers to isolate sections of chemical scheme based on MCM KPP format
    drh_str = params['drh_str']
    erh_str = params['erh_str']
    const_comp = params['const_comp']
    Init_comp = params['Init_comp']
    Diff_setname = params['Diff_setname']
    Diff_set = params['Diff_set']
    con_infl_nam = const_comp
    Zgrid = params['Zgrid']  # number of grid points in tube length direction
    Rgrid = params['Rgrid']  # number of grid points in tube radius direction
    # formula = params['formula'] # the formula for the plots
    key_spe_for_plot = params['key_spe_for_plot']  # key species for ploting
    plot_spec = params['plot_spec']  # plot species
    dt = params['dt']
    ratio = params['ratio']
    H2Oconc = const_comp_conc[:, const_comp.index('H2O')]

    # read the file and store everything into a list
    f_open_eqn = open(sch_name, mode='r')  # open the chemical scheme file
    total_list_eqn = f_open_eqn.readlines()
    f_open_eqn.close()  # close file

    eqn_list, num_eqn, rrc, rrc_name, RO2_names, eqn_list_on = sch_interr.sch_interr(total_list_eqn, chm_sch_mrk)

    [rindx, rstoi, pindx, pstoi, reac_coef,
     nreac, nprod, y_arr, y_rind, uni_y_rind, y_pind,
     uni_y_pind, reac_col, prod_col, rstoi_flat, pstoi_flat,
     rr_arr, rr_arr_p, comp_namelist, comp_list,
     comp_num] = eqn_interr.eqn_interr(num_eqn, eqn_list, chm_sch_mrk)

    [rindx, pindx, rstoi, pstoi, nreac, nprod, y_arr, y_rind,
     uni_y_rind, y_pind, uni_y_pind, reac_col, prod_col,
     rstoi_flat, pstoi_flat, rr_arr, rr_arr_p,
     comp_num, RO2_indx,
     HOMRO2_indx, comp_list,
     eqn_num, comp_namelist,
     erf, err_mess, con_C_indx] = eqn_pars.extr_mech(sch_name, chm_sch_mrk,
                                                     con_infl_nam, const_comp,
                                                     drh_str, erh_str)

    RO2_indi = RO2_indices.RO2_indices(comp_namelist, RO2_names)
    u, Diff_vals = get_diff_and_u(comp_namelist, Diff_setname, con_C_indx, Diff_set, T, p)
    numLoop = 500  # number of times to run to reach the pinhole of the instrument
    timesteps = 1000  # number of timesteps, dt * timesteps * numLoop is time elapsed in the final solution
    # print (comp_namelist, Diff_vals)
    # Change odd number Rgrid to even number grid
    if (Rgrid % 2) != 0:
        Rgrid = Rgrid + 1

    sp_line = int(Zgrid * L1 / (L2 + L1))
    # % set the dr dx Q parameters for the tube
    Qtot = np.zeros([int(Rgrid), int(Zgrid), comp_num])
    Qtot[:, 0:sp_line, :] = Q1
    Qtot[:, sp_line:, :] = Q2

    Rtot = np.zeros([int(Rgrid), int(Zgrid), comp_num])
    Rtot[:, 0:sp_line, :] = R1
    Rtot[:, sp_line:, :] = R2

    const_comp_gird = cal_const_comp_conc.cal_const_comp_conc(Rgrid, Zgrid, const_comp_conc, L1, L2, const_comp)
    const_comp_grid_1 = cal_const_comp_conc.cal_const_comp_conc_1(Rgrid, Zgrid, const_comp_conc, L1, L2, const_comp)
    c = np.zeros([Rgrid, Zgrid, comp_num])
    for i in Init_comp:
        c[:, 0, comp_namelist.index(i)] = Init_comp_conc[
            Init_comp.index(i)]  # set [OH] at z = 0 # set [HO2] at z = 0. This equals OH conc

    write_rate_file.write_rate_file(reac_coef, p, rrc, rrc_name, 0)

    C0 = c[0, 0, :]

    [y, y_mw, num_comp, M, y_indx_plot, dydt_vst,
     comp_namelist, erf, err_mess] = init_conc.init_conc(comp_num, comp_namelist, C0, T, \
                                                         p, comp_namelist, rindx, pindx, \
                                                         num_eqn[0], nreac, nprod, comp_namelist, \
                                                         RO2_indx, HOMRO2_indx, rstoi, pstoi)
    # y = c[0,0,:]

    import rate_coeffs
    RO2conc = RO2_conc.RO2_conc(RO2_indi, y)
    op = jude_species(y, comp_namelist)

    rate_values, erf, err_mess = rate_coeffs.evaluate_rates(RO2conc, T, 0, M, M * 0.7809, op[0], op[1], op[2], op[3],
                                                            op[4], p)

    # % plot
    if Q2 > Q1:
        # first tube run
        for i in const_comp:
            c[:, :, comp_namelist.index(i)] = const_comp_grid_1[const_comp.index(i)]
        dr = np.zeros([int(Rgrid), int(Zgrid), comp_num])
        dx = L1 / (Zgrid - 1)
        dr[:, :, :] = 2 * R1 / (Rgrid - 1)
        for j in range(numLoop):
            c1 = c.copy()
            old = c1[:, -1, comp_namelist.index(key_spe_for_plot)]

            c = odesolve_single(timesteps, Zgrid, Rgrid, dt, Diff_vals, R1, dr, dx, Q1, c, comp_namelist, dydt_vst, rindx,
                         nreac, rstoi, rate_values, const_comp, u)

            new = c[:, -1, comp_namelist.index(key_spe_for_plot)]

            if (j > 5) & (np.sum(new - old) / np.sum(old) < 1e-5):
                break
        # second tube run
        dr = np.zeros([int(Rgrid), int(Zgrid), comp_num])
        dx = L2 / (Zgrid - 1)
        dr[:, :, :] = 2 * R2 / (Rgrid - 1)
        for j in range(numLoop):
            c1 = c.copy()
            old = c1[:, -1, comp_namelist.index(key_spe_for_plot)]

            c = odesolve_Y(timesteps, Zgrid, Rgrid, dt, Diff_vals, R2, dr, dx, Q2, c, comp_namelist, dydt_vst, rindx,
                           nreac, rstoi, rate_values, const_comp, u)

            tim = (j + 1) * timesteps * dt
            comp_plot_index = [comp_namelist.index(plot_spec[i]) for i in range(len(plot_spec))]

            new = c[:, -1, comp_namelist.index(key_spe_for_plot)]

            fig, axs = plt.subplots(2, 3, figsize=(8, 5), facecolor='w', edgecolor='k')
            fig.subplots_adjust(hspace=.5, wspace=.45)
            plt.style.use('default')
            plt.rcParams.update(
                {'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif', 'font.serif': 'Times New Roman'})

            axs = axs.ravel()
            formula = get_formula(plot_spec)
            for i in range(len(plot_spec)):
                #axs[i].pcolor(np.linspace(0, L2 + L1, Zgrid), np.linspace(-R1, R1, Rgrid), c[:, :, comp_plot_index[i]],
                 #             shading='nearest', cmap='jet')
                axs[i].pcolor(np.linspace(0, L2, Zgrid), np.linspace(-R2, R2, Rgrid), c[:, :, comp_plot_index[i]],
                              shading='nearest', cmap='jet')

                axs[i].set_xlabel('L [cm]')
                axs[i].set_ylabel('R [cm]')
                axs[i].set_title(formula[i])

            fig.delaxes(axs[5])
            plt.gcf().text(0.7, 0.3, 'Time = ' + str(tim), fontsize=15)
            plt.draw()
            plt.pause(1)

            print(['t = ' + str(tim) + str(key_spe_for_plot) + " difference: " + str(np.sum(new - old))])

            if (j > 5) & (np.sum(new - old) / np.sum(old) < 1e-5):
                break

        dr_final = R2 / (Rgrid - 1) * 2
        x = np.arange(0, R2, dr_final) + dr_final
        rVec = np.arange(0, R2, 0.001)

        meanConc = []
        for i in plot_spec:
            y_x = np.flip(c[0: int(Rgrid / 2), -1, comp_namelist.index(i)])  # 'SA'
            splineres1 = interpolate.splrep(x, y_x)
            cVec = interpolate.splev(rVec, splineres1)
            meanConc.append(2 * 0.001 / R2 ** 2 * np.sum(cVec * rVec))
    else:
        for i in const_comp:
            c[:, :, comp_namelist.index(i)] = const_comp_gird[const_comp.index(i)]
        if R2 == 0:
            R2 = R1
        dr = np.zeros([int(Rgrid), int(Zgrid), comp_num])
        dx = (L2 + L1) / (Zgrid - 1)
        dr[:, 0:sp_line, :] = 2 * R1 / (Rgrid - 1)
        dr[:, sp_line:, :] = 2 * R2 / (Rgrid - 1)
        for j in range(numLoop):
            c1 = c.copy()
            old = c1[:, -1, comp_namelist.index(key_spe_for_plot)]

            c = odesolve_double(timesteps, Zgrid, Rgrid, dt, Diff_vals, Rtot, dr, dx, Qtot, c, comp_namelist, dydt_vst,
                           rindx,
                           nreac, rstoi, rate_values, const_comp, u)

            tim = (j + 1) * timesteps * dt
            comp_plot_index = [comp_namelist.index(plot_spec[i]) for i in range(len(plot_spec))]

            new = c[:, -1, comp_namelist.index(key_spe_for_plot)]

            fig, axs = plt.subplots(2, 3, figsize=(8, 5), facecolor='w', edgecolor='k')
            fig.subplots_adjust(hspace=.5, wspace=.45)
            plt.style.use('default')
            plt.rcParams.update(
                {'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif', 'font.serif': 'Times New Roman'})

            axs = axs.ravel()
            formula = get_formula(plot_spec)
            for i in range(len(plot_spec)):
                axs[i].pcolor(np.linspace(0, L2 + L1, Zgrid), np.linspace(-R1, R1, Rgrid), c[:, :, comp_plot_index[i]],
                              shading='nearest', cmap='jet')
                axs[i].pcolor(np.linspace(0, L2 + L1, Zgrid), np.linspace(-R1, R2, Rgrid), c[:, :, comp_plot_index[i]],
                              shading='nearest', cmap='jet')

                axs[i].set_xlabel('L [cm]')
                axs[i].set_ylabel('R [cm]')
                axs[i].set_title(formula[i])

            fig.delaxes(axs[5])
            plt.gcf().text(0.7, 0.3, 'Time = ' + str(tim), fontsize=15)
            plt.draw()
            plt.pause(1)

            print(['t = ' + str(tim) + str(key_spe_for_plot) + " difference: " + str(np.sum(new - old))])

            if (j > 5) & (np.sum(new - old) / np.sum(old) < 1e-5):
                break
        dr_final = R2 / (Rgrid - 1) * 2
        x = np.arange(0, R2, dr_final) + dr_final
        rVec = np.arange(0, R2, 0.001)

        meanConc = []
        for i in plot_spec:
            y_x = np.flip(c[0: int(Rgrid / 2), -1, comp_namelist.index(i)])  # 'SA'

            splineres1 = interpolate.splrep(x, y_x)

            cVec = interpolate.splev(rVec, splineres1)
            meanConc.append(2 * 0.001 / R2 ** 2 * np.sum(cVec * rVec))

    return meanConc, c
