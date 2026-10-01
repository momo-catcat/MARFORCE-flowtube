import numpy as np
import pandas as pd
import multiprocessing
import os
import signal
from Funcs.meanconc_cal import meanconc_cal_sim as meanconc_cal_sim
from Funcs.meanconc_cal import meanconc_cal
from Funcs.odesolve3 import odesolve as odesolve
from Funcs.ode_solv_batch import ode_solv_batch
from Funcs.cal_const_comp_conc import set_const_comp_conc_for_1st_tube, set_const_comp_conc_for_2nd_tube
from Funcs.set_Rdrdx_sepera_simu import set_Rtot_dr_dx_for_spe_simulation
from Funcs.plots_updating import plot_concentration_profiles, plot_concentration_profiles_3
from Funcs.set_boundlay_for_ode import (set_boundary_conditions,
                                        set_boundary_conditions_point,
                                        set_boundary_conditions_continuousOH_point)
from Funcs.ode_worker import init_worker, build_initargs
import csv

# Checkpoint save interval in seconds (30 min)
_CKPT_INTERVAL = 1800


def _ckpt_path(modelparams):
    """Build checkpoint file path — unique per parameter combination."""
    return (f"{modelparams.export_file_folder}"
            f"R{modelparams.Rgrid}L{modelparams.Zgrid}"
            f"inter_acc{modelparams.actual_acc:.0e}"
            f"dt{modelparams.dt:.0e}timstep{modelparams.timesteps:.0e}"
            f"itx{modelparams.Itx:.1e}{modelparams.model_mode}"
            f"_Stage_{modelparams.number_stage}_checkpoint.npz")


def _save_checkpoint(path, c, c_2nd, k, delta_c_final, tim_1_final, final_sp,
                     comp_namelist):
    """Save iteration state to disk so the run can be resumed (no-op when path is None)."""
    if path is None:
        return
    # Convert final_sp dict-of-lists to a 2D array for efficient storage
    sp_arr = np.array([final_sp[spe] for spe in comp_namelist], dtype=np.float64)
    np.savez_compressed(path,
                        c=c, c_2nd=c_2nd,
                        k=np.int64(k),
                        delta_c_final=np.array(delta_c_final),
                        tim_1_final=np.array(tim_1_final),
                        sp_arr=sp_arr)
    print(f'  [checkpoint] saved iter {k} -> {os.path.basename(path)}', flush=True)


def _load_checkpoint(path, comp_namelist):
    """Load checkpoint. Returns (c, c_2nd, k_start, delta_c_final, tim_1_final, final_sp) or None."""
    if path is None or not os.path.isfile(path):
        return None
    try:
        data = np.load(path, allow_pickle=False)
        c = data['c']
        c_2nd = data['c_2nd']
        k_start = int(data['k']) + 1  # resume from next iteration
        delta_c_final = data['delta_c_final'].tolist()
        tim_1_final = data['tim_1_final'].tolist()
        sp_arr = data['sp_arr']
        final_sp = {spe: sp_arr[i].tolist() for i, spe in enumerate(comp_namelist)}
        print(f'  [checkpoint] loaded iter {k_start - 1} from {os.path.basename(path)}', flush=True)
        return c, c_2nd, k_start, delta_c_final, tim_1_final, final_sp
    except Exception as e:
        print(f'  [checkpoint] failed to load ({e}); starting fresh', flush=True)
        return None


def model_twotubes(numLoop, Diff_vals, rowvals, colptrs, u, plot_spec, formula, c, Q1, Q2, modelparams):
    comp_plot_indices = [modelparams.comp_namelist.index(plot_spec[i]) for i in range(len(plot_spec))]
    key_index = modelparams.comp_namelist.index(modelparams.key_spe_for_plot)
    old_1 = np.zeros([modelparams.Rgrid, modelparams.Zgrid, modelparams.comp_num])

    # Integration step for chemistry (constant across iterations)
    integ_step = float(modelparams.dt) * int(modelparams.timesteps)

    c_2nd = np.zeros([modelparams.Rgrid2, modelparams.Zgrid2, modelparams.comp_num])
    # ---- Warm-start c_2nd from previous run or previous stage ----
    _ws_c2 = getattr(modelparams, '_ws_c2nd', None)
    if _ws_c2 is not None and _ws_c2.shape == c_2nd.shape:
        c_2nd[:] = _ws_c2
        print(f'  [warm start] loaded c_2nd ({c_2nd.shape}) from previous state', flush=True)
    if hasattr(modelparams, '_ws_c2nd'):
        del modelparams._ws_c2nd
    dr1 = modelparams.R1 / (modelparams.Rgrid1 - 1)
    dr2 = modelparams.R2 / (modelparams.Rgrid2 - 1)
    delta_c_final = []
    tim_1_final = []
    # Species mean-concentration snapshots recorded every iteration (matches reference)
    final_sp = {spe: [] for spe in modelparams.comp_namelist}

    # Precompute second-tube spatial parameters (constant across iterations)
    dr_2nd, dx_2nd = set_Rtot_dr_dx_for_spe_simulation(modelparams.R2, modelparams.L2, modelparams)
    OHsource_2nd = 'point'

    # Precompute diffusion arrays outside the loop (they never change)
    D = np.broadcast_to(Diff_vals, (modelparams.Rgrid, modelparams.Zgrid, modelparams.comp_num))
    D2 = np.broadcast_to(Diff_vals, (modelparams.Rgrid2, modelparams.Zgrid2, modelparams.comp_num))

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
    ckpt_file = _ckpt_path(modelparams) if modelparams.use_restart_now else None
    k_start = 0
    ckpt = _load_checkpoint(ckpt_file, modelparams.comp_namelist)
    if ckpt is not None:
        c, c_2nd, k_start, delta_c_final, tim_1_final, final_sp = ckpt

    # ---- SIGTERM handler: save checkpoint before SLURM kills us ----
    _sigterm_received = [False]

    def _sigterm_handler(signum, frame):
        _sigterm_received[0] = True
        print(f'\n  [checkpoint] SIGTERM received — will save and exit after current iter',
              flush=True)

    if ckpt_file is not None:
        signal.signal(signal.SIGTERM, _sigterm_handler)

    _last_ckpt_time = _time.perf_counter()

    for k in range(k_start, numLoop):
        _t0 = _time.perf_counter()
        tim = round((k + 1) * modelparams.timesteps * modelparams.dt, 3)

        # ---- First tube -------------------------------------------------------
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
                c_a = c_reshaped[indices_a]
                c_b = c_reshaped[indices_b]

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
            c = (set_boundary_conditions_continuousOH_point(c, old_1, modelparams)
                 if modelparams.Init_set == 'on'
                 else set_boundary_conditions(c, old_1, modelparams))
            c = odesolve(modelparams.timesteps, modelparams.Zgrid, modelparams.Rgrid,
                         modelparams.dt, D, modelparams.R1, modelparams.dr, modelparams.dx,
                         Q1, c, u, modelparams.rrc_calculated, modelparams, modelparams.OHsource)

        else:
            c = odesolve(modelparams.timesteps, modelparams.Zgrid, modelparams.Rgrid,
                         modelparams.dt, D, modelparams.R1, modelparams.dr, modelparams.dx,
                         Q1, c, u, modelparams.rrc_calculated, modelparams, modelparams.OHsource)

        # ---- Second tube ------------------------------------------------------
        for i in range(modelparams.Rgrid2):
            c_2nd[i, 0, :] = (c[i, -1, :] * Q1 / Q2
                               * ((2 * dr1 * i * dr1 + dr1 ** 2) / (modelparams.R1 ** 2))
                               / ((2 * dr2 * i * dr2 + dr2 ** 2) / modelparams.R2 ** 2))
        c_2nd = set_const_comp_conc_for_2nd_tube(c_2nd, modelparams)

        old = c_2nd[:, -1, key_index].copy()
        old_1_2nd = c_2nd.copy()

        if 'flowtube2' in modelparams.model_mode:
            c2_reshaped = c_2nd.reshape(-1, modelparams.comp_num)
            c2_reshaped = ode_solv_batch(c2_reshaped, integ_step,
                                         modelparams.rrc_calculated_p,
                                         modelparams, rowvals, colptrs)
            c_2nd = c2_reshaped.reshape(modelparams.Rgrid2, modelparams.Zgrid2, modelparams.comp_num)

            c_2nd = set_const_comp_conc_for_2nd_tube(c_2nd, modelparams)
            c_2nd = set_boundary_conditions_point(c_2nd, old_1_2nd, modelparams)
            c_2nd = odesolve(modelparams.timesteps, modelparams.Zgrid2, modelparams.Rgrid2,
                             modelparams.dt, D2, modelparams.R2, dr_2nd, dx_2nd,
                             Q2, c_2nd, u, modelparams.rrc_calculated_p, modelparams, OHsource_2nd)
        else:
            c_2nd = odesolve(modelparams.timesteps, modelparams.Zgrid2, modelparams.Rgrid2,
                             modelparams.dt, D2, modelparams.R2, dr_2nd, dx_2nd,
                             Q2, c_2nd, u, modelparams.rrc_calculated_p, modelparams, OHsource_2nd)

        new = c_2nd[:, -1, key_index]
        old_norm = np.linalg.norm(old)
        if old_norm == 0:
            delta_c = np.inf
        else:
            delta_c = np.linalg.norm(new - old) / old_norm
        tim_1 = round((k + 1) * modelparams.timesteps * modelparams.dt, 3)
        delta_c_final.append(delta_c)
        tim_1_final.append(tim_1)

        # Tube 1 delta: how much tube 1 exit changed this iteration
        _t1_new = c[:, -1, key_index]
        _t1_old_norm = np.linalg.norm(old_1[:, -1, key_index])
        _t1_delta = np.linalg.norm(_t1_new - old_1[:, -1, key_index]) / _t1_old_norm if _t1_old_norm > 0 else float('inf')

        _elapsed = _time.perf_counter() - _t0
        print(f'  iter {k}: t={tim_1:.3f}s  delta={delta_c:.6f}  t1_delta={_t1_delta:.6f}  elapsed={_elapsed:.2f}s', flush=True)

        # NaN guard: abort only on true NaN (never on Inf — first iter sets delta_c=inf
        # legitimately when old_norm==0, so checking !isfinite would abort valid runs).
        if np.isnan(delta_c) or np.isnan(c_2nd).any():
            print(f'  iter {k}: NaN detected — aborting this run '
                  f'(likely CFL-unstable: dt={modelparams.dt:.0e} too large for this grid)',
                  flush=True)
            delta_c_final.pop()
            tim_1_final.pop()
            break

        # Record mean concentrations every iteration — one call for all species
        all_means = meanconc_cal(c_2nd, modelparams)
        for spe, val in zip(modelparams.comp_namelist, all_means):
            final_sp[spe].append(val)

        if k % modelparams.num_plot == 0:
            if modelparams.OHsource == 'point':
                plot_concentration_profiles(c, plot_spec, formula,
                                            (modelparams.L2 + modelparams.L1), 0,
                                            modelparams.Zgrid, modelparams.Rgrid, tim,
                                            modelparams.R2, comp_plot_indices)
            else:
                plot_concentration_profiles_3(c, c_2nd, plot_spec, formula, tim,
                                              comp_plot_indices, modelparams)
            print('t = ' + str(tim) + ' ' + str(modelparams.key_spe_for_plot) + " difference ratio: {:.4f} ".format(delta_c)+str(modelparams.key_spe_for_plot)+" conc: {:.2E}".format(new[int(modelparams.Rgrid2 / 2)]))

        # ---- Periodic checkpoint save (every _CKPT_INTERVAL seconds) ----
        _now = _time.perf_counter()
        if _now - _last_ckpt_time >= _CKPT_INTERVAL:
            _save_checkpoint(ckpt_file, c, c_2nd, k, delta_c_final,
                             tim_1_final, final_sp, modelparams.comp_namelist)
            _last_ckpt_time = _now

        # ---- SIGTERM: save checkpoint and exit cleanly ----
        if _sigterm_received[0]:
            _save_checkpoint(ckpt_file, c, c_2nd, k, delta_c_final,
                             tim_1_final, final_sp, modelparams.comp_namelist)
            print(f'  [checkpoint] exiting after iter {k} due to SIGTERM', flush=True)
            break

        if k > modelparams.fix_timstep and delta_c < _conv_threshold:
            # Converged — remove checkpoint file (no longer needed)
            if ckpt_file is not None and os.path.isfile(ckpt_file):
                os.remove(ckpt_file)
                print(f'  [checkpoint] converged — removed {os.path.basename(ckpt_file)}',
                      flush=True)
            break

    # Restore default SIGTERM handler
    signal.signal(signal.SIGTERM, signal.SIG_DFL)

    delta_base = (f"{modelparams.export_file_folder}"
                  f"R{modelparams.Rgrid}L{modelparams.Zgrid}"
                  f"inter_acc{modelparams.actual_acc:.0e}"
                  f"dt{modelparams.dt:.0e}timstep{modelparams.timesteps:.0e}"
                  f"itx{modelparams.Itx:.1e}{modelparams.model_mode}"
                  f"_Stage_{modelparams.number_stage}")

    # Save ONE combined CSV: time + delta + all species mean concs (matches reference)
    df_dict = {"time": tim_1_final, "delta": delta_c_final}
    df_dict.update(final_sp)
    pd.DataFrame(df_dict).to_csv(delta_base + "_delta_keyspecforplot.csv", index=False)

    output_file = (f"{modelparams.export_file_folder}"
                   f"R{modelparams.Rgrid}L{modelparams.Zgrid}"
                   f"inter_acc{modelparams.actual_acc:.0e}"
                   f"dt{modelparams.dt:.0e}timstep{modelparams.timesteps:.0e}"
                   f"itx{modelparams.Itx:.1e}{modelparams.model_mode}"
                   f"_Stage_{modelparams.number_stage}_1st_tube_profile_output.txt")
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        write = csv.writer(f)
        write.writerows(c)

    # Cleanup multiprocessing pool
    if getattr(modelparams, '_pool', None) is not None:
        modelparams._pool.close()
        modelparams._pool.join()
        modelparams._pool = None

    return c, c_2nd
