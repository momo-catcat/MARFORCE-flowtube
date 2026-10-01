
def set_boundary_conditions(c, old_1, modelparams):
    # Vectorized: use precomputed index array of variable (non-const) components.
    # modelparams.var_comp_indices is set once in cmd_calib5 after eqn_pars.extr_mech.
    idx = modelparams.var_comp_indices
    c[0, :, idx] = old_1[0, :, idx]
    c[-1, :, idx] = old_1[-1, :, idx]
    return c


def set_boundary_conditions_continuousOH_point(c, old_1, modelparams):
    idx = modelparams.var_comp_indices
    c[0, :, idx] = old_1[0, :, idx]
    c[-1, :, idx] = old_1[-1, :, idx]

    # Restore z=0 inlet values for Init_comp species (non-const only).
    # modelparams.var_init_comp_indices is precomputed in cmd_calib5.
    init_idx = modelparams.var_init_comp_indices
    if len(init_idx):
        c[1:-1, 0, init_idx] = old_1[1:-1, 0, init_idx]
    return c


def set_boundary_conditions_point(c, old_1, modelparams):
    idx = modelparams.var_comp_indices
    c[0, :, idx] = old_1[0, :, idx]
    c[-1, :, idx] = old_1[-1, :, idx]
    c[1:-1, 0, idx] = old_1[1:-1, 0, idx]
    return c
