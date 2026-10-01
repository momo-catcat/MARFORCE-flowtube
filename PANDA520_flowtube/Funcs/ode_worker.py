"""Pool-worker module for multiprocessing-based cell-level ODE solving.

The key optimisation here is the *initializer* pattern:
  - All large, constant arrays (mechanism indices, stoichiometries, Jacobian
    structure, etc.) are loaded into each worker process **once** at pool
    creation time via init_worker().
  - Per-iteration, only the small per-cell concentration vector ``y`` and the
    reaction-rate-coefficient array ``rrc`` need to be sent over IPC.
  - This eliminates the repeated pickling of hundreds of kB of constant data
    that would otherwise happen on every pool.starmap() call.
"""

import Funcs.ode_solv as _ode_solv_mod

# Module-level dict filled once by the pool initializer.
_C = {}


def init_worker(rindx_g, pindx_g, rstoi_g, pstoi_g, nreac_g, nprod_g,
                jac_stoi_g, njac_g, jac_den_indx_g, jac_indx_g,
                y_arr_g, y_rind_g, uni_y_rind_g, y_pind_g, uni_y_pind_g,
                reac_col_g, prod_col_g, rstoi_flat_g, pstoi_flat_g,
                rr_arr_g, rr_arr_p_g, rowvals, colptrs, comp_num,
                dil_fac_now, wall_loss, integ_step,
                jac_flat_rr_idx, jac_flat_stoi,
                jac_flat_den_indx, jac_flat_data_indx):
    """Called once per worker process to pre-load all constant data."""
    global _C
    _C = dict(
        rindx_g=rindx_g, pindx_g=pindx_g, rstoi_g=rstoi_g, pstoi_g=pstoi_g,
        nreac_g=nreac_g, nprod_g=nprod_g,
        jac_stoi_g=jac_stoi_g, njac_g=njac_g,
        jac_den_indx_g=jac_den_indx_g, jac_indx_g=jac_indx_g,
        y_arr_g=y_arr_g, y_rind_g=y_rind_g, uni_y_rind_g=uni_y_rind_g,
        y_pind_g=y_pind_g, uni_y_pind_g=uni_y_pind_g,
        reac_col_g=reac_col_g, prod_col_g=prod_col_g,
        rstoi_flat_g=rstoi_flat_g, pstoi_flat_g=pstoi_flat_g,
        rr_arr_g=rr_arr_g, rr_arr_p_g=rr_arr_p_g,
        rowvals=rowvals, colptrs=colptrs, comp_num=comp_num,
        dil_fac_now=dil_fac_now, wall_loss=wall_loss, integ_step=integ_step,
        jac_flat_rr_idx=jac_flat_rr_idx, jac_flat_stoi=jac_flat_stoi,
        jac_flat_den_indx=jac_flat_den_indx, jac_flat_data_indx=jac_flat_data_indx,
    )


def solve_cell(args):
    """Solve the chemistry ODE for a single grid cell.

    Parameters
    ----------
    args : tuple (y, rrc)
        y   – 1-D concentration array for this cell (molecules/cm³)
        rrc – reaction-rate coefficients for this section/stage
              (sent per call because it can differ between sections)
    """
    y, rrc = args
    c = _C
    return _ode_solv_mod.ode_solv(
        y, c['integ_step'], rrc,
        c['rindx_g'], c['pindx_g'], c['rstoi_g'], c['pstoi_g'],
        c['nreac_g'], c['nprod_g'],
        c['jac_stoi_g'], c['njac_g'], c['jac_den_indx_g'], c['jac_indx_g'],
        c['y_arr_g'], c['y_rind_g'], c['uni_y_rind_g'],
        c['y_pind_g'], c['uni_y_pind_g'],
        c['reac_col_g'], c['prod_col_g'], c['rstoi_flat_g'], c['pstoi_flat_g'],
        c['rr_arr_g'], c['rr_arr_p_g'],
        c['rowvals'], c['colptrs'], c['comp_num'],
        c['dil_fac_now'], c['wall_loss'],
        c['jac_flat_rr_idx'], c['jac_flat_stoi'],
        c['jac_flat_den_indx'], c['jac_flat_data_indx'],
    )


def build_initargs(modelparams, rowvals, colptrs):
    """Helper – collect all constant arrays into the tuple expected by init_worker."""
    return (
        modelparams.rindx_g, modelparams.pindx_g,
        modelparams.rstoi_g, modelparams.pstoi_g,
        modelparams.nreac_g, modelparams.nprod_g,
        modelparams.jac_stoi_g, modelparams.njac_g,
        modelparams.jac_den_indx_g, modelparams.jac_indx_g,
        modelparams.y_arr_g, modelparams.y_rind_g, modelparams.uni_y_rind_g,
        modelparams.y_pind_g, modelparams.uni_y_pind_g,
        modelparams.reac_col_g, modelparams.prod_col_g,
        modelparams.rstoi_flat_g, modelparams.pstoi_flat_g,
        modelparams.rr_arr_g, modelparams.rr_arr_p_g,
        rowvals, colptrs, modelparams.comp_num,
        modelparams.dil_fac_now, modelparams.wall_loss,
        float(modelparams.dt) * int(modelparams.timesteps),
        modelparams.jac_flat_rr_idx, modelparams.jac_flat_stoi,
        modelparams.jac_flat_den_indx, modelparams.jac_flat_data_indx,
    )
