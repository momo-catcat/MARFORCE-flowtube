"""Numbalsoda + numba.prange parallel chemistry ODE batch solver.

Design
------
All mechanism arrays (indices, stoichiometries, rate coefficients) are packed
into one flat float64 array `_packed_data`.  A Numba @cfunc function
`_lsoda_rhs` reads from this pointer — no Python callbacks, no GIL.

`numba_batch_solve` uses `numba.prange` (OpenMP-style threads) to solve all
grid cells in parallel.  No inter-process IPC; threads share the read-only
packed mechanism data and each write to their own slice of Y_new.

Packed-data layout
------------------
  [0]  n_react  (int cast to float64)
  [1]  n_max
  [2]  n_comp
  [3]  n_pow    (len of pow_mask_flat)
  [4]  n_loss   (len of y_rind_g)
  [5]  n_prod   (len of y_pind_g)
  [6]  dil_fac_now
  [7 .. 7+n_comp-1]                         wall_loss[n_comp]
  [7+n_comp .. +n_react-1]                  rrc[n_react]        ← updated each iter
  [+n_loss]                                 y_arr_g[n_loss]
  [+n_loss]                                 y_rind_g[n_loss]
  [+n_loss]                                 rr_arr_g[n_loss]
  [+n_loss]                                 rstoi_flat_g[n_loss]
  [+n_prod]                                 y_pind_g[n_prod]
  [+n_prod]                                 rr_arr_p_g[n_prod]
  [+n_prod]                                 pstoi_flat_g[n_prod]
  [+n_pow]                                  pow_mask_flat[n_pow]
  [+n_pow]                                  pow_stoi[n_pow]
"""

import numpy as np
import numba
from numba import njit, cfunc, prange
from numba import carray as _carray
from numbalsoda import lsoda, lsoda_sig, address_as_void_pointer


# ── Packed-data helpers ────────────────────────────────────────────────────────

def pack_mechanism_data(rrc, y_arr_g, y_rind_g, rr_arr_g, rstoi_flat_g,
                        y_pind_g, rr_arr_p_g, pstoi_flat_g,
                        pow_mask_flat, pow_stoi,
                        dil_fac_now, wall_loss,
                        n_react, n_max, n_comp):
    """Build and return the packed float64 data array (once per run).

    Returns (packed_data, off_rrc) where `off_rrc` is the offset of the rrc
    section so it can be updated cheaply each outer iteration.
    """
    n_pow  = len(pow_mask_flat)
    n_loss = len(y_rind_g)
    n_prod = len(y_pind_g)

    header = np.array([n_react, n_max, n_comp, n_pow, n_loss, n_prod,
                       float(dil_fac_now)], dtype=np.float64)

    packed = np.concatenate([
        header,                                  # 7
        np.asarray(wall_loss,     dtype=np.float64),   # n_comp
        np.asarray(rrc,           dtype=np.float64)[:n_react],  # n_react
        np.asarray(y_arr_g,       dtype=np.float64),   # n_loss
        np.asarray(y_rind_g,      dtype=np.float64),   # n_loss
        np.asarray(rr_arr_g,      dtype=np.float64),   # n_loss
        np.asarray(rstoi_flat_g,  dtype=np.float64),   # n_loss
        np.asarray(y_pind_g,      dtype=np.float64),   # n_prod
        np.asarray(rr_arr_p_g,    dtype=np.float64),   # n_prod
        np.asarray(pstoi_flat_g,  dtype=np.float64),   # n_prod
        np.asarray(pow_mask_flat, dtype=np.float64),   # n_pow
        np.asarray(pow_stoi,      dtype=np.float64),   # n_pow
    ])
    off_rrc = int(7 + n_comp)
    return packed, off_rrc


def update_rrc_in_packed(packed_data, rrc, off_rrc, n_react):
    """Update the rrc section of packed_data in-place (cheap, called each iter)."""
    packed_data[off_rrc: off_rrc + n_react] = rrc[:n_react]


# ── Numba @cfunc RHS (no Python, no GIL) ─────────────────────────────────────

@cfunc(lsoda_sig)
def _lsoda_rhs(t, u, du, p):
    """LSODA RHS: reads all mechanism data from the packed pointer p."""
    # Header
    n_react = int(p[0])
    n_max   = int(p[1])
    n_comp  = int(p[2])
    n_pow   = int(p[3])
    n_loss  = int(p[4])
    n_prod  = int(p[5])
    dil_fac = p[6]

    # Offsets (same order as pack_mechanism_data)
    off_wall   = 7
    off_rrc    = off_wall  + n_comp
    off_yarr   = off_rrc   + n_react
    off_yrind  = off_yarr  + n_loss
    off_rrarr  = off_yrind + n_loss
    off_rstoi  = off_rrarr + n_loss
    off_ypind  = off_rstoi + n_loss
    off_rrp    = off_ypind + n_prod
    off_pstoi  = off_rrp   + n_prod
    off_pmask  = off_pstoi + n_prod
    off_pstoi2 = off_pmask + n_pow

    # Build rrc_y (flat, shape n_react×n_max), initialised to 1
    rrc_y = np.ones(n_react * n_max)
    for k in range(n_loss):
        y_arr_k  = int(p[off_yarr  + k])
        y_rind_k = int(p[off_yrind + k])
        rrc_y[y_arr_k] = u[y_rind_k]

    # Power corrections
    for k in range(n_pow):
        idx = int(p[off_pmask + k])
        rrc_y[idx] = rrc_y[idx] ** p[off_pstoi2 + k]

    # Reaction rates
    rr = np.empty(n_react)
    for i in range(n_react):
        val = p[off_rrc + i]
        for j in range(n_max):
            val *= rrc_y[i * n_max + j]
        rr[i] = val

    # Initialise du
    for i in range(n_comp):
        du[i] = 0.0

    # Loss of reactants
    for k in range(n_loss):
        rr_arr_k  = int(p[off_rrarr + k])
        y_rind_k  = int(p[off_yrind + k])
        rstoi_k   = p[off_rstoi  + k]
        du[y_rind_k] -= rr[rr_arr_k] * rstoi_k

    # Gain of products
    for k in range(n_prod):
        rr_p_k   = int(p[off_rrp   + k])
        y_pind_k = int(p[off_ypind + k])
        pstoi_k  = p[off_pstoi + k]
        du[y_pind_k] += rr[rr_p_k] * pstoi_k

    # Dilution + wall loss
    for i in range(n_comp):
        du[i] -= u[i] * dil_fac
        du[i] -= u[i] * p[off_wall + i]


# ── Parallel batch solver ──────────────────────────────────────────────────────

@njit(parallel=True, cache=False)   # can't cache: ctypes pointer inside lsoda
def numba_batch_solve(Y, integ_step, packed_data, funcptr,
                      rtol=1e-4, atol=1e-5):
    """Solve chemistry ODE for all N_cells using numba.prange (thread parallel).

    Parameters
    ----------
    Y           : float64[N_cells, comp_num]   initial concentrations
    integ_step  : float                         integration step (s)
    packed_data : float64[DATA_SIZE]            packed mechanism (read-only)
    funcptr     : int                           address of _lsoda_rhs @cfunc
    rtol, atol  : float                         LSODA tolerances

    Returns
    -------
    Y_new : float64[N_cells, comp_num]
    """
    N_cells = Y.shape[0]
    comp_num = Y.shape[1]
    Y_new = np.empty_like(Y)
    t_eval = np.array([0.0, integ_step])

    for i in prange(N_cells):
        u0 = Y[i].copy()
        usol, success = lsoda(address_as_void_pointer(funcptr),
                               u0, t_eval, data=packed_data,
                               rtol=rtol, atol=atol)
        if success:
            Y_new[i] = usol[1]   # solution at t = integ_step
        else:
            # Flag failure the same way ode_solv does (first element = -1e6)
            Y_new[i] = u0
            Y_new[i, 0] = -1.0e6
    return Y_new


# ── Module-level state (set once per run) ────────────────────────────────────

_packed_data = None   # set by init_numba_batch()
_off_rrc     = None
_n_react     = None
_funcptr     = int(_lsoda_rhs.address)   # C function pointer address


def init_numba_batch(modelparams, rrc, rowvals, colptrs):
    """Initialise packed data for the current mechanism (call once per stage)."""
    global _packed_data, _off_rrc, _n_react

    mp = modelparams
    n_react = int(mp.rindx_g.shape[0])
    n_max   = int(mp.rindx_g.shape[1])
    n_comp  = int(mp.comp_num)

    # Power-correction flat indices
    pow_mask2d   = (mp.rstoi_g != 1.0) & (mp.rstoi_g != 0.0)
    pow_mask_flat = np.flatnonzero(pow_mask2d).astype(np.int64)
    pow_stoi      = mp.rstoi_g.ravel()[pow_mask_flat]

    _packed_data, _off_rrc = pack_mechanism_data(
        rrc,
        mp.y_arr_g, mp.y_rind_g, mp.rr_arr_g, mp.rstoi_flat_g,
        mp.y_pind_g, mp.rr_arr_p_g, mp.pstoi_flat_g,
        pow_mask_flat, pow_stoi,
        mp.dil_fac_now, mp.wall_loss,
        n_react, n_max, n_comp,
    )
    _n_react = n_react
    print(f'  [numba_batch] packed data: {len(_packed_data)} float64 values '
          f'({len(_packed_data)*8/1024:.1f} KB)', flush=True)


def solve_batch_numba(Y, integ_step, rrc):
    """Top-level entry called from ode_solv_batch.py.

    Updates rrc in packed data, then runs parallel LSODA.
    """
    global _packed_data, _off_rrc, _n_react, _funcptr
    if _packed_data is None:
        raise RuntimeError('init_numba_batch() must be called before solve_batch_numba()')

    # Update reaction rate coefficients (the only part that changes each iter)
    update_rrc_in_packed(_packed_data, rrc, _off_rrc, _n_react)

    return numba_batch_solve(Y, integ_step, _packed_data, _funcptr)
