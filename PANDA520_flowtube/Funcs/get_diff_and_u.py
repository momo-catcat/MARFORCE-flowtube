import numpy as np
from kinetics.diff_coef import diff_coef
from Funcs.diffusion_const_added import add_diff_const as add_diff_const

def get_diff_and_u(comp_namelist, Diff_setname, con_C_indx, Diff_set, T, p, modelparams):
    # Compute diffusion coefficients
    Diff_vals = np.array([diff_coef(i, T, p, 'air', modelparams) for i in comp_namelist], dtype=np.float32)

    # Replace specific species with predefined diffusion values
    for i in Diff_setname:
        Diff_vals[comp_namelist.index(i)] = np.float32(Diff_set[Diff_setname.index(i)])

    # Get indices of non-constant compounds
    u = np.array([i for j, i in enumerate(range(len(comp_namelist))) if j not in con_C_indx], dtype=np.int32)

    return u, Diff_vals


## This function is used because of the O1D here
def get_diff_and_u_for_more_species(comp_namelist, Diff_setname, con_C_indx, Diff_set, T, p, modelparams):
    modelparams.Diff_vals = np.full(len(comp_namelist), np.nan, dtype=np.float32)  # Initialize as NaNs

    for i in comp_namelist:
        idx = comp_namelist.index(i)
        if i == 'H2O':
            modelparams.Diff_vals[idx] = np.float32(add_diff_const(p, T)[1])
        elif i not in Diff_setname:
            modelparams.Diff_vals[idx] = np.float32(diff_coef(i, T, p, 'air', modelparams)[0])
        else:
            modelparams.Diff_vals[idx] = np.float32(Diff_set[Diff_setname.index(i)])

    # Get indices of non-constant compounds
    modelparams.u = np.array([i for j, i in enumerate(range(len(comp_namelist))) if j not in con_C_indx], dtype=np.int32)

    return modelparams.u, modelparams.Diff_vals