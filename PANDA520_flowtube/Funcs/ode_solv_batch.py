"""Batch chemistry ODE solver: solves ALL grid cells in parallel.

Priority order
--------------
1. Numba+numbalsoda path  — numba.prange threads, no GIL, no IPC overhead.
   Enabled after init_numba_batch() has been called (done in cmd_calib5).
2. Multiprocessing pool   — separate processes, good for non-Numba envs.
3. Sequential fallback    — single-threaded, used when pool creation failed.
"""

import numpy as np
from Funcs.ode_solv import ode_solv
from Funcs.ode_worker import solve_cell
from Funcs.ode_solv_numba_batch import solve_batch_numba
import Funcs.ode_solv_numba_batch as _nb_mod

_pool_checked   = False  # print path status once on first call
_numba_checked  = False


def ode_solv_batch(Y, integ_step, rrc, modelparams, rowvals, colptrs):
    """Solve chemistry ODE for ALL grid cells by looping over cells.

    Parameters
    ----------
    Y : ndarray, shape (N_cells, comp_num)
        Initial concentrations for all grid cells (molecules/cm³).
    integ_step : float
        Integration time step (s).
    rrc : ndarray
        Reaction rate coefficients array.
    modelparams : object
        Must carry all mechanism arrays needed by ode_solv.
        If modelparams._pool exists, uses multiprocessing.
    rowvals : ndarray
        Row indices of Jacobian non-zero elements.
    colptrs : ndarray
        Column pointers for the sparse Jacobian.

    Returns
    -------
    Y_new : ndarray, shape (N_cells, comp_num)
        Updated concentrations after integ_step seconds.
    """
    global _pool_checked, _numba_checked
    N_cells  = Y.shape[0]
    comp_num = int(modelparams.comp_num)

    # ---------- Numba+numbalsoda path (highest priority) ----------
    if _nb_mod._packed_data is not None:
        if not _numba_checked:
            print(f'  [ode_solv_batch] NUMBA parallel path active '
                  f'({N_cells} cells, numba.prange)', flush=True)
            _numba_checked = True
        try:
            return solve_batch_numba(Y.astype(np.float64), integ_step, rrc)
        except Exception as e:
            print(f'  [ode_solv_batch] numba path failed ({e}), falling back',
                  flush=True)
            _nb_mod._packed_data = None   # disable for remaining iterations

    # ---------- Multiprocessing path ---------
    pool = getattr(modelparams, '_pool', None)
    if not _pool_checked:
        if pool is not None:
            print(f'  [ode_solv_batch] PARALLEL pool path active ({N_cells} cells)', flush=True)
        else:
            print(f'  [ode_solv_batch] SEQUENTIAL path active (no pool)', flush=True)
        _pool_checked = True
    if pool is not None:
        try:
            args = [(Y[i].astype(float), rrc) for i in range(N_cells)]
            results = pool.map(solve_cell, args, chunksize=max(1, N_cells // (4 * 10)))
            return np.array(results, dtype=float)
        except Exception as e:
            print(f'  [ode_solv_batch] pool.map failed ({e}), falling back to sequential', flush=True)
            modelparams._pool = None

    # ---------- Sequential fallback ----------
    Y_new    = np.empty((N_cells, comp_num), dtype=float)

    # Unpack all mechanism arrays once (avoid repeated attribute lookups)
    rindx_g          = modelparams.rindx_g
    pindx_g          = modelparams.pindx_g
    rstoi_g          = modelparams.rstoi_g
    pstoi_g          = modelparams.pstoi_g
    nreac_g          = modelparams.nreac_g
    nprod_g          = modelparams.nprod_g
    jac_stoi_g       = modelparams.jac_stoi_g
    njac_g           = modelparams.njac_g
    jac_den_indx_g   = modelparams.jac_den_indx_g
    jac_indx_g       = modelparams.jac_indx_g
    y_arr_g          = modelparams.y_arr_g
    y_rind_g         = modelparams.y_rind_g
    uni_y_rind_g     = modelparams.uni_y_rind_g
    y_pind_g         = modelparams.y_pind_g
    uni_y_pind_g     = modelparams.uni_y_pind_g
    reac_col_g       = modelparams.reac_col_g
    prod_col_g       = modelparams.prod_col_g
    rstoi_flat_g     = modelparams.rstoi_flat_g
    pstoi_flat_g     = modelparams.pstoi_flat_g
    rr_arr_g         = modelparams.rr_arr_g
    rr_arr_p_g       = modelparams.rr_arr_p_g
    dil_fac_now      = modelparams.dil_fac_now
    wall_loss        = np.asarray(modelparams.wall_loss, dtype=float)
    jac_flat_rr_idx  = modelparams.jac_flat_rr_idx
    jac_flat_stoi    = modelparams.jac_flat_stoi
    jac_flat_den_indx = modelparams.jac_flat_den_indx
    jac_flat_data_indx = modelparams.jac_flat_data_indx

    for i in range(N_cells):
        Y_new[i] = ode_solv(
            Y[i].astype(float), integ_step, rrc,
            rindx_g, pindx_g, rstoi_g, pstoi_g, nreac_g, nprod_g,
            jac_stoi_g, njac_g, jac_den_indx_g, jac_indx_g,
            y_arr_g, y_rind_g, uni_y_rind_g, y_pind_g, uni_y_pind_g,
            reac_col_g, prod_col_g, rstoi_flat_g, pstoi_flat_g, rr_arr_g,
            rr_arr_p_g, rowvals, colptrs, comp_num, dil_fac_now, wall_loss,
            jac_flat_rr_idx, jac_flat_stoi, jac_flat_den_indx, jac_flat_data_indx
        )

    return Y_new

