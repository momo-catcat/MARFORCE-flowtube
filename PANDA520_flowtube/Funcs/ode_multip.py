import Funcs.ode_solv

def solve_ode(y, dt, rrc, rindx_g, pindx_g, rstoi_g, pstoi_g, nreac_g, nprod_g,
              jac_stoi_g, njac_g, jac_den_indx_g, jac_indx_g,
              y_arr_g, y_rind_g, uni_y_rind_g, y_pind_g, uni_y_pind_g,
              reac_col_g, prod_col_g, rstoi_flat_g, pstoi_flat_g, rr_arr_g,
              rr_arr_p_g, rowvals, colptrs, comp_num, dil_fac_now=0, wall_loss=0,
              jac_flat_rr_idx=None, jac_flat_stoi=None,
              jac_flat_den_indx=None, jac_flat_data_indx=None):
    y = Funcs.ode_solv.ode_solv(
        y, dt, rrc,
        rindx_g, pindx_g, rstoi_g, pstoi_g, nreac_g, nprod_g,
        jac_stoi_g, njac_g, jac_den_indx_g, jac_indx_g,
        y_arr_g, y_rind_g, uni_y_rind_g, y_pind_g, uni_y_pind_g,
        reac_col_g, prod_col_g, rstoi_flat_g, pstoi_flat_g, rr_arr_g,
        rr_arr_p_g, rowvals, colptrs, comp_num, dil_fac_now, wall_loss,
        jac_flat_rr_idx, jac_flat_stoi, jac_flat_den_indx, jac_flat_data_indx)
    return y
