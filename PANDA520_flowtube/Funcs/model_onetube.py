import numpy as np
import pandas as pd
import multiprocessing
import os
import signal
from Funcs.meanconc_cal import meanconc_cal_sim as meanconc_cal_sim
from Funcs.meanconc_cal import meanconc_cal
from Funcs.odesolve3 import odesolve as odesolve
from Funcs.ode_solv_batch import ode_solv_batch
from Funcs.cal_const_comp_conc import set_const_comp_conc_for_1st_tube
from Funcs.plots_updating import plot_concentration_profiles, plot_concentration_profiles_different_dx
from Funcs.set_boundlay_for_ode import set_boundary_conditions, set_boundary_conditions_point
from Funcs.ode_worker import init_worker, build_initargs
from Funcs.model_twotubes import _ckpt_path, _save_checkpoint, _load_checkpoint, _CKPT_INTERVAL


def model_onetube(numLoop, Diff_vals, rowvals, colptrs, u, plot_spec, formula, c, Q1, Q2, modelparams):
    comp_plot_indices = [modelparams.comp_namelist.index(plot_spec[i]) for i in range(len(plot_spec))]
    key_index = modelparams.comp_namelist.index(modelparams.key_spe_for_plot)
    old_1 = np.zeros([modelparams.Rgrid, modelparams.Zgrid, modelparams.comp_num])

    # Integration step for chemistry (same value that was passed to Pool workers)
    integ_step = float(modelparams.dt) * int(modelparams.timesteps)

    delta_c_final = []
    tim_1_final = []
    # Species mean-concentration snapshots recorded every iteration (matches reference)
    final_sp = {spe: [] for spe in modelparams.comp_namelist}

    # Create multiprocessing pool for parallel cell ODE solves
    n_workers = os.cpu_count() or 4
    try:
        _pool = multiprocessing.Pool(
            processes=n_workers,
            initializer=init_worker,
            initargs=build_initargs(modelparams, rowvals, colptrs),
        )
        modelparams._pool = _pool
        print(f'  Pool created with {n_workers} workers', flush=True)
    except Exception as e:
        print(f'  Pool creation failed ({e}), running sequentially', flush=True)
        modelparams._pool = None

    import time as _time

    # ---- Convergence threshold: allow override via modelparams.converge_threshold ----
    _conv_threshold = getattr(modelparams, 'converge_threshold',
                              modelparams.actual_acc_1st)
    if _conv_threshold != modelparams.actual_acc_1st:
        print(f'  [converge] using override threshold {_conv_threshold:.0e} '
              f'(actual_acc_1st={modelparams.actual_acc_1st:.0e})', flush=True)

    # ---- Checkpoint: load previous state if available ----
    ckpt_file = _ckpt_path(modelparams)
    k_start = 0
    # model_onetube has no c_2nd; use a dummy for the shared checkpoint API
    _dummy_c2 = np.empty(0)
    ckpt = _load_checkpoint(ckpt_file, modelparams.comp_namelist)
    if ckpt is not None:
        c, _dummy_c2, k_start, delta_c_final, tim_1_final, final_sp = ckpt

    # ---- SIGTERM handler ----
    _sigterm_received = [False]

    def _sigterm_handler(signum, frame):
        _sigterm_received[0] = True
        print(f'\n  [checkpoint] SIGTERM received — will save and exit after current iter',
              flush=True)

    signal.signal(signal.SIGTERM, _sigterm_handler)
    _last_ckpt_time = _time.perf_counter()

    for k in range(k_start, numLoop):
        tim = round((k + 1) * modelparams.timesteps * modelparams.dt, 3)
        old = c[:, -1, key_index].copy()
        D = np.broadcast_to(Diff_vals, (modelparams.Rgrid, modelparams.Zgrid, modelparams.comp_num))

        if 'flowtube2' in modelparams.model_mode:
            old_1 = c.copy()

            if modelparams.OHsource == 'Continuous':
                indices_a = [i * modelparams.Zgrid + j
                             for i in range(modelparams.Rgridl)
                             for j in range(modelparams.Zgridl)]
                indices_b = [i * modelparams.Zgrid + j
                             for i in range(modelparams.Rgrid1)
                             for j in range(modelparams.Zgridl, modelparams.Zgrid)]

                c_reshaped = c.reshape(-1, modelparams.comp_num)
                c_a = c_reshaped[indices_a]   # (Rgridl*Zgridl, comp_num)
                c_b = c_reshaped[indices_b]   # (Rgrid1*Zgrid1, comp_num)

                # Batch solve each section with its own rrc
                c_a = ode_solv_batch(c_a, integ_step, modelparams.rrc_calculated,
                                     modelparams, rowvals, colptrs)
                c_a = c_a.reshape(modelparams.Rgridl, modelparams.Zgridl, modelparams.comp_num)

                c_b = ode_solv_batch(c_b, integ_step, modelparams.rrc_calculated_p,
                                     modelparams, rowvals, colptrs)
                c_b = c_b.reshape(modelparams.Rgrid1, modelparams.Zgrid1, modelparams.comp_num)

                c = np.concatenate((c_a, c_b), axis=1)

            else:
                c_reshaped = c.reshape(-1, modelparams.comp_num)
                c_reshaped = ode_solv_batch(c_reshaped, integ_step,
                                            modelparams.rrc_calculated,
                                            modelparams, rowvals, colptrs)
                c = c_reshaped.reshape(modelparams.Rgrid, modelparams.Zgrid, modelparams.comp_num)

            c = set_const_comp_conc_for_1st_tube(c, modelparams)
            c = (set_boundary_conditions_point(c, old_1, modelparams)
                 if modelparams.Init_set == 'on'
                 else set_boundary_conditions(c, old_1, modelparams))
            c = odesolve(modelparams.timesteps, modelparams.Zgrid, modelparams.Rgrid,
                         modelparams.dt, D, modelparams.R1, modelparams.dr, modelparams.dx,
                         Q1, c, u, modelparams.rrc_calculated, modelparams, modelparams.OHsource)

        else:
            c = odesolve(modelparams.timesteps, modelparams.Zgrid, modelparams.Rgrid,
                         modelparams.dt, D, modelparams.R1, modelparams.dr, modelparams.dx,
                         Q1, c, u, modelparams.rrc_calculated, modelparams, modelparams.OHsource)

        new = c[:, -1, key_index]
        old_norm = np.linalg.norm(old)
        if old_norm == 0:
            delta_c = np.inf
        else:
            delta_c = np.linalg.norm(new - old) / old_norm
        tim_1 = round((k + 1) * modelparams.timesteps * modelparams.dt, 3)
        delta_c_final.append(delta_c)
        tim_1_final.append(tim_1)

        # NaN guard: abort only on true NaN (never on Inf — first iter sets delta_c=inf
        # legitimately when old_norm==0, so checking !isfinite would abort valid runs).
        if np.isnan(delta_c) or np.isnan(c).any():
            print(f'  iter {k}: NaN detected — aborting this run '
                  f'(likely CFL-unstable: dt={modelparams.dt:.0e} too large for this grid)',
                  flush=True)
            break

        # Record mean concentrations every iteration — one call for all species
        all_means = meanconc_cal(c, modelparams)
        for spe, val in zip(modelparams.comp_namelist, all_means):
            final_sp[spe].append(val)

        if k % modelparams.num_plot == 0:
            if modelparams.OHsource == 'point':
                plot_concentration_profiles(c, plot_spec, formula,
                                            (modelparams.L2 + modelparams.L1), 0,
                                            modelparams.Zgrid, modelparams.Rgrid, tim,
                                            modelparams.R2, comp_plot_indices)
            else:
                plot_concentration_profiles_different_dx(c, plot_spec, formula, tim,
                                                         comp_plot_indices, modelparams)
            print('t = ' + str(tim) + ' ' + str(modelparams.key_spe_for_plot) + " difference ratio: {:.4f} ".format(delta_c)+str(modelparams.key_spe_for_plot)+" conc: {:.2E}".format(new[int(modelparams.Rgrid / 2)]))

        # ---- Periodic checkpoint save ----
        _now = _time.perf_counter()
        if _now - _last_ckpt_time >= _CKPT_INTERVAL:
            _save_checkpoint(ckpt_file, c, _dummy_c2, k, delta_c_final,
                             tim_1_final, final_sp, modelparams.comp_namelist)
            _last_ckpt_time = _now

        # ---- SIGTERM: save checkpoint and exit ----
        if _sigterm_received[0]:
            _save_checkpoint(ckpt_file, c, _dummy_c2, k, delta_c_final,
                             tim_1_final, final_sp, modelparams.comp_namelist)
            print(f'  [checkpoint] exiting after iter {k} due to SIGTERM', flush=True)
            break

        if k > modelparams.fix_timstep and delta_c < _conv_threshold:
            if os.path.isfile(ckpt_file):
                os.remove(ckpt_file)
                print(f'  [checkpoint] converged — removed {os.path.basename(ckpt_file)}',
                      flush=True)
            break

    signal.signal(signal.SIGTERM, signal.SIG_DFL)

    # Cleanup multiprocessing pool
    if getattr(modelparams, '_pool', None) is not None:
        modelparams._pool.close()
        modelparams._pool.join()
        modelparams._pool = None

    # Save ONE combined CSV: time + delta + all species mean concs (matches reference)
    delta_path = (f"{modelparams.export_file_folder}"
                  f"R{modelparams.Rgrid}L{modelparams.Zgrid}"
                  f"inter_acc{modelparams.actual_acc:.0e}"
                  f"dt{modelparams.dt:.0e}timstep{modelparams.timesteps:.0e}"
                  f"itx{modelparams.Itx:.1e}{modelparams.model_mode}"
                  f"_Stage_{modelparams.number_stage}_delta_keyspecforplot.csv")

    df_dict = {"time": tim_1_final, "delta": delta_c_final}
    df_dict.update(final_sp)
    df = pd.DataFrame(df_dict)
    df.to_csv(delta_path, index=False)

    return c
