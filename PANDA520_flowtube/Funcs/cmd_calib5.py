import numpy as np
import os
import Funcs.eqn_pars as eqn_pars
import Funcs.init_conc as init_conc
import Funcs.RO2_indices as RO2_indices
from Funcs.cal_const_comp_conc import set_const_comp_conc_for_1st_tube
from Funcs.judg_spe_reac_rates import jude_species as jude_species
from Funcs.get_diff_and_u import get_diff_and_u_for_more_species
from Funcs.get_formula import get_formula
from Funcs.model_onetube import model_onetube
from Funcs.model_twotubes import model_twotubes
from Funcs.model_box import model_box
from Funcs.meanconc_cal import meanconc_cal
from Funcs.grid_parameters import grid_para as grid_para
import Funcs.group_indices as group_indices
import Funcs.rrc_calc as rrc_calc
from Funcs.wall_loss_cal_species import get_wall_loss
from Funcs.ode_solv_numba_batch import init_numba_batch

def cmd_calib5( const_comp_conc, modelparams,Init_comp_conc, Q1, Q2, c_prev=None):
    modelparams.stage_const_comp_conc = const_comp_conc
    Q1 = Q1 / 60  # flow for first tube
    Q2 = Q2  / 60  # flow for second tube
    modelparams.chem_sch_mrk = ['<', 'RO2', '+', 'C(ind_', ')','' , '&', '' , '', ':', '>', ';', '','%']
    modelparams.light_stat_now = 0 
    modelparams.con_infl_nam = modelparams.const_comp
    modelparams.comp0 = modelparams.Init_comp
    # Restart files (checkpoints, warm start, stage carry-over) are meant for long continuous-OH runs.
    # Point-source calibrations always start fresh; override with use_restart = True/False in the start script.
    use_restart = bool(getattr(modelparams, 'use_restart', modelparams.OHsource == 'Continuous'))
    # flow-tube transport is solved in the flowtube modes and in 'kinetic' mode (transport only, no chemistry)
    tube_mode = ('flowtube' in modelparams.model_mode) or (modelparams.model_mode == 'kinetic')
    modelparams.use_restart_now = use_restart
    # print('Q1', modelparams.Q1, 'Q2', modelparams.Q2)
    try: ### if the O2 concentration is not given, set it to 0
        modelparams.O2 = modelparams.O2conc[0][0]
    except:
        modelparams.O2 = 0

    # read the file and separate the equations and rate coefficients
    int_tol = 1e-5,1e-4
    [rrc, rrc_name, rowvals, colptrs, comp_num, Jlen, erf, err_mess, modelparams] = eqn_pars.extr_mech(int_tol, 0, modelparams)

    # ------------------------------------------------------------------
    # Precompute flat Jacobian arrays for vectorised jac() in ode_solv.
    # These are constant for the entire run (fixed mechanism), so we build
    # them once here and store in modelparams instead of rebuilding inside
    # every solve_ivp call.
    _nreac = int(modelparams.eqn_num[0])
    _njac  = modelparams.njac_g[:_nreac, 0].astype(int)
    _row   = np.repeat(np.arange(_nreac), _njac)              # reaction index per flat element
    _col   = np.concatenate([np.arange(n) for n in _njac])    # local column index per element
    modelparams.jac_flat_rr_idx    = _row
    modelparams.jac_flat_stoi      = modelparams.jac_stoi_g[:_nreac][_row, _col]
    modelparams.jac_flat_den_indx  = modelparams.jac_den_indx_g[:_nreac][_row, _col].astype(int)
    modelparams.jac_flat_data_indx = modelparams.jac_indx_g[:_nreac][_row, _col].astype(int)

    # Precompute variable-component index array for vectorised boundary conditions.
    _const_set = set(modelparams.const_comp)
    modelparams.var_comp_indices = np.array(
        [i for i, name in enumerate(modelparams.comp_namelist) if name not in _const_set],
        dtype=int)
    # Indices of Init_comp species that are NOT constant-concentration
    modelparams.var_init_comp_indices = np.array(
        [modelparams.comp_namelist.index(s)
         for s in modelparams.Init_comp if s not in _const_set],
        dtype=int)

    # Precompute net stoichiometry matrix for batch chemistry solver.
    # rxn_net_stoich[i, s] = net change in species s per unit rate of reaction i.
    # Shape: (nreac, comp_num). Built once here; used in ode_solv_batch.
    _rxn_net = np.zeros((_nreac, comp_num), dtype=np.float64)
    np.add.at(_rxn_net,
              (modelparams.rr_arr_p_g, modelparams.y_pind_g),
              modelparams.pstoi_flat_g)
    np.add.at(_rxn_net,
              (modelparams.rr_arr_g, modelparams.y_rind_g),
              -modelparams.rstoi_flat_g)
    modelparams.rxn_net_stoich = _rxn_net
    # ------------------------------------------------------------------

    modelparams = RO2_indices.RO2_indices(modelparams)

    for i in range(len(modelparams.RO2_indices)):
        nu = modelparams.RO2_indices[i][1]
        # print(modelparams.Pybel_objects[nu].formula)
        modelparams = group_indices.group_indices(modelparams.comp_smil, nu, modelparams)
    
    # get the diffusion for all species and  the index of species in C except constant compounds
    u, Diff_vals = get_diff_and_u_for_more_species(modelparams.comp_namelist, modelparams.Diff_setname, modelparams.con_C_indx, modelparams.Diff_set, modelparams.TEMP, modelparams.p, modelparams)
    numLoop = 500000000  # number of times to run to reach the pinhole of the instrument
    if modelparams.runtime[modelparams.number_stage] /(modelparams.dt * modelparams.timesteps) > modelparams.fix_timstep:
        modelparams.fix_timstep = modelparams.runtime[modelparams.number_stage] /(modelparams.dt * modelparams.timesteps)
        modelparams.actual_acc_1st = 1
    else:
        modelparams.actual_acc_1st = modelparams.inter_acc

    # Change odd number Rgrid to even number grid
    if (modelparams.Rgrid % 2) != 0:
        modelparams.Rgrid = modelparams.Rgrid + 1
    # Two-tube model when the flow changes (Y-piece) or a 2nd section with a different radius follows
    modelparams.two_tubes = bool(Q1 != Q2 or (float(getattr(modelparams, 'L2', 0)) > 0
                                               and float(modelparams.R2) != float(modelparams.R1)))
    modelparams = grid_para(modelparams)
    # print(modelparams.comp_namelist)
    # % set the concentration for all the species for the grid of 80*40 in c
    c = np.zeros([modelparams.Rgrid, modelparams.Zgrid, comp_num])


    y = np.zeros([comp_num])
    for i in modelparams.const_comp:  # set the constant concentrations for const_comp
        try:
            y[modelparams.comp_namelist.index(i)] = modelparams.stage_const_comp_conc[0, modelparams.const_comp.index(i)]
            c[:,:,modelparams.comp_namelist.index(i)] = modelparams.stage_const_comp_conc[0, modelparams.const_comp.index(i)]
        except:
            pass
    if modelparams.Init_set == 'on':
        for i in modelparams.comp0:  # set [OH] at z = 0 # set [HO2] at z = 0 [oh]
            c[:, 0, modelparams.comp_namelist.index(i)] = Init_comp_conc[modelparams.comp0.index(i)]
            y[modelparams.comp_namelist.index(i)] = Init_comp_conc[modelparams.comp0.index(i)]
    c_box = y
    (y, H2Oi, y_mw, num_comp, Cfactor, y_indx_plot, corei, inj_indx,  nrec_steps, erf, err_mess, NOi, HO2i, NO3i, init_conc_c, modelparams) = init_conc.init_conc(modelparams.H2Oconc[0][0], comp_num, c[0, 0, :], modelparams.p, 0, modelparams.eqn_num[0], modelparams.comp0, modelparams)

    # used as the title for the plotted figures
    modelparams = get_formula(modelparams)
    
    # % set the grids parameters

    [rrc_calculated, erf, err_mess] = rrc_calc.rrc_calc(modelparams.H2Oconc[modelparams.number_stage][0], modelparams.TEMP, y, modelparams.p, Jlen, y[NOi], y[HO2i], y[NO3i], 0, modelparams)

    for index, reaction in enumerate(modelparams.eqn_list):
        if "itProd" in reaction:
            modelparams.itprod_index = index
            modelparams.rrc_calculated = rrc_calculated
            modelparams.rrc_calculated_p = rrc_calculated.copy()
            modelparams.rrc_calculated_p[modelparams.itprod_index] = 0
        else:
            modelparams.rrc_calculated = rrc_calculated
            modelparams.rrc_calculated_p = rrc_calculated.copy()  

    ######### set constant concentration for the whole grid 
    c = set_const_comp_conc_for_1st_tube(c,modelparams)

    ######### prepare dilution rate for box model
    if hasattr(modelparams, 'chamber_volume') and hasattr(modelparams, 'total_flow'):
        # chamber_volume in liters, total_flow in slpm → k_dil in s⁻¹
        modelparams.dil_fac_now = np.float32(
            (modelparams.total_flow / 60.0) / modelparams.chamber_volume)
        print(f'  [box] dilution rate = {modelparams.dil_fac_now:.4e} s⁻¹ '
              f'(V={modelparams.chamber_volume} L, Q={modelparams.total_flow} slpm)')

    ######### prepare wall loss for box model
    if modelparams.wall_loss_set == 1:
        modelparams.wall_loss = get_wall_loss(modelparams.TEMP,modelparams.p, modelparams)
    else:
        modelparams.wall_loss = np.ones([comp_num], dtype=np.float32) * modelparams.wall_loss_set

    # Apply user-defined wall loss for specific species
    if hasattr(modelparams, 'wall_loss_custom'):
        for spe, kw in modelparams.wall_loss_custom.items():
            if spe in modelparams.comp_namelist:
                modelparams.wall_loss[modelparams.comp_namelist.index(spe)] = kw

    modelparams.wall_loss = np.array(modelparams.wall_loss, dtype=float)
    if np.isinf(modelparams.wall_loss).any():
        modelparams.wall_loss[np.isinf(modelparams.wall_loss)] = 2e-3

    # ---- Initialise Numba+numbalsoda parallel batch solver (#9) ----
    if 'flowtube' in modelparams.model_mode:
        try:
            init_numba_batch(modelparams, modelparams.rrc_calculated,
                             rowvals, colptrs)
        except Exception as _e:
            print(f'  [cmd_calib5] Numba batch init failed ({_e}); '
                  f'falling back to pool/sequential', flush=True)

    # ---- Warm start (#6): load converged state from a previous run if available ----
    _ws_loaded = False
    _ws_base = (f"{modelparams.export_file_folder}"
                f"warmstart_stage{modelparams.number_stage}_"
                f"R{modelparams.Rgrid}L{modelparams.Zgrid}")
    _ws_npz = _ws_base + '.npz'
    _ws_npy = _ws_base + '.npy'

    if use_restart and os.path.isfile(_ws_npz):
        # New format: contains both first tube (c) and second tube (c_2nd)
        _ws = np.load(_ws_npz)
        if 'c' in _ws and _ws['c'].shape == c.shape:
            c = _ws['c']
            _ws_loaded = True
            print(f'  [warm start] loaded c {c.shape} from {_ws_npz}', flush=True)
        if 'c_2nd' in _ws:
            modelparams._ws_c2nd = _ws['c_2nd']
            print(f'  [warm start] loaded c_2nd {_ws["c_2nd"].shape} from {_ws_npz}', flush=True)
    elif use_restart and os.path.isfile(_ws_npy):
        # Old format: single array (c_2nd only, shape mismatch with c)
        _c_ws = np.load(_ws_npy)
        if _c_ws.shape == c.shape:
            c = _c_ws
            _ws_loaded = True
            print(f'  [warm start] loaded {_ws_npy}', flush=True)
        elif (_c_ws.shape[0] == c.shape[0] and _c_ws.shape[2] == c.shape[2]
              and _c_ws.shape[1] < c.shape[1]):
            # Old saved data is c_2nd with fewer z-points.
            # Map into main-tube section of c, and also pass to model_twotubes.
            _zl = getattr(modelparams, 'Zgridl', 0)
            _zs = _c_ws.shape[1]
            if _zl + _zs <= c.shape[1]:
                c[:, _zl:_zl + _zs, :] = _c_ws
                modelparams._ws_c2nd = _c_ws
                _ws_loaded = True
                print(f'  [warm start] partial load {_c_ws.shape} → '
                      f'c[:, {_zl}:{_zl + _zs}, :] + c_2nd from {_ws_npy}', flush=True)
            else:
                print(f'  [warm start] shape mismatch {_c_ws.shape} vs {c.shape}; '
                      f'ignoring', flush=True)
        else:
            print(f'  [warm start] shape mismatch {_c_ws.shape} vs {c.shape}; '
                  f'ignoring', flush=True)

    # ---- Stage carry-over: use previous stage's converged state if no warm-start ----
    if use_restart and not _ws_loaded and c_prev is not None:
        if isinstance(c_prev, list):
            # c_prev is a list of 1D arrays from box model — use last converged state
            c_box = c_prev[-1].copy()
            print(f'  [stage init] loaded previous box-model state ({len(c_prev)} iterations)', flush=True)
        elif isinstance(c_prev, tuple):
            # c_prev = (c_1st, c_2nd) from twotubes
            _c1, _c2 = c_prev
            if _c1.shape == c.shape:
                c[:] = _c1
                print(f'  [stage init] loaded previous stage c {_c1.shape}', flush=True)
            modelparams._ws_c2nd = _c2
            print(f'  [stage init] loaded previous stage c_2nd {_c2.shape}', flush=True)
        elif c_prev.shape == c.shape:
            c[:] = c_prev
            print(f'  [stage init] loaded previous stage converged state', flush=True)
        elif (c_prev.shape[0] == c.shape[0] and c_prev.shape[2] == c.shape[2]
              and c_prev.shape[1] < c.shape[1]):
            _zl = getattr(modelparams, 'Zgridl', 0)
            _zs = c_prev.shape[1]
            if _zl + _zs <= c.shape[1]:
                c[:, _zl:_zl + _zs, :] = c_prev
                modelparams._ws_c2nd = c_prev
                print(f'  [stage init] loaded previous stage → '
                      f'c[:, {_zl}:{_zl + _zs}, :] + c_2nd', flush=True)

    # A loaded field carries the previous stage's inlet values: re-apply this stage's inlet
    if tube_mode and modelparams.Init_set == 'on':
        for i in modelparams.comp0:
            c[:, 0, modelparams.comp_namelist.index(i)] = Init_comp_conc[modelparams.comp0.index(i)]

    ########################################################################################### run the model
    tube_messages = {
        1: 'one diameter tube',
        2: 'two diameter tubes',
        3: 'two diameter tubes, with different flow rates'
    }

    flag = int(modelparams.flag_tube)  # Ensure it's an int
    print(f'tube_flag {flag}, {tube_messages.get(flag, "unknown configuration")}')

    print(
        'OH source is continuous, not point'
        if modelparams.OHsource == 'Continuous'
        else 'OH source is point, not continuous'
    )

    if tube_mode:
        if modelparams.two_tubes:
            c, c_2nd = model_twotubes(numLoop, Diff_vals, rowvals, colptrs,u, modelparams.plot_spec, modelparams.formula, c, Q1,Q2,modelparams)
            meanConc = meanconc_cal(c_2nd, modelparams)          # second tube → uses R2 (default)
        else:
            c = model_onetube(numLoop, Diff_vals, rowvals, colptrs,u, modelparams.plot_spec, modelparams.formula, c, Q1,Q2,modelparams)
            meanConc = meanconc_cal(c, modelparams, R=modelparams.R1)  # single tube → use R1
  
    elif 'box' in modelparams.model_mode:
        print('box model')
        # Zero wall_loss and dilution for constant species — they are held
        # constant by resetting after each integration step, but the ODE solver
        # must not apply loss terms to them *during* integration, otherwise
        # they decay and all downstream chemistry sees artificially low
        # concentrations of the held species.
        _const_set = set(modelparams.const_comp)
        for idx_c, name_c in enumerate(modelparams.comp_namelist):
            if name_c in _const_set:
                modelparams.wall_loss[idx_c] = 0.0
        # Convert scalar dilution rate to per-species array so constant
        # species get zero dilution inside the ODE solver.
        _dil_arr = np.full(comp_num, float(modelparams.dil_fac_now))
        for idx_c, name_c in enumerate(modelparams.comp_namelist):
            if name_c in _const_set:
                _dil_arr[idx_c] = 0.0
        modelparams.dil_fac_now = _dil_arr

        c = model_box(numLoop,  modelparams.key_spe_for_plot, modelparams.dt, rowvals, colptrs,modelparams.plot_spec, modelparams.formula, c_box, rrc_calculated,modelparams)
        meanConc = c[-1]

    # % print the meanconc for key_spe_for_plot
    print( modelparams.key_spe_for_plot + "meanconc: {:.2E}".format(meanConc[modelparams.comp_namelist.index(modelparams.key_spe_for_plot)]))
    print( "OH meanconc: {:.2E}".format(meanConc[modelparams.comp_namelist.index('OH')]))

    # ---- Warm start (#6): save converged state for future runs ----
    if tube_mode:
        _ws_base = (f"{modelparams.export_file_folder}"
                    f"warmstart_stage{modelparams.number_stage}_"
                    f"R{modelparams.Rgrid}L{modelparams.Zgrid}")
        if modelparams.two_tubes:
            # Two-tube: save both first tube (c) and second tube (c_2nd)
            if use_restart:
                np.savez(_ws_base + '.npz', c=c, c_2nd=c_2nd)
                print(f'  [warm start] saved c {c.shape} + c_2nd {c_2nd.shape} '
                      f'→ {_ws_base}.npz', flush=True)
            return meanConc, (c, c_2nd)
        else:
            # One-tube: save c directly
            if use_restart:
                np.save(_ws_base + '.npy', c)
                print(f'  [warm start] saved {_ws_base}.npy', flush=True)
            return meanConc, c

    return meanConc, c
