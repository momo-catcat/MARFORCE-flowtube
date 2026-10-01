# %import packages
import numpy as np
from numba import njit


@njit(cache=True)
def _odesolve_inner(timesteps, Zgrid, Rgrid, dt, D, Rtot, dr, dx,
                    Qtot, c, u, sp):
    """Numba-accelerated inner loop for diffusion + convection PDE.

    Matches the reference vectorized algorithm exactly:
    - Half-tube (rows 0..sp-1) computed, then mirrored to sp..Rgrid-1
    - Row 0 left unchanged (term1/term2 are zero there)
    - Non-uniform grid second derivative for axial diffusion
    - 1/r radial term: (-1/r) * (c[ir]-c[ir-1]) / dr
    """
    dr2 = dr * dr
    conv = (2.0 * Qtot) / (np.pi * Rtot ** 4)
    num_u = len(u)

    # Precompute radial positions (absolute distance from centre)
    r1d = np.empty(Rgrid, dtype=np.float64)
    for ir in range(Rgrid):
        r1d[ir] = abs((ir - Rgrid / 2.0 + 0.5) * dr)

    # Pre-allocate snapshot buffer once — avoids one heap allocation per timestep
    initc = np.empty_like(c)

    for m in range(timesteps):
        # In-place copy into pre-allocated buffer (no allocation)
        initc[:] = c

        # --- interior points (radial 1..sp-1, axial 1..Zgrid-2) ---
        for ir in range(1, sp):
            rv = r1d[ir]
            flow_fac = conv * (Rtot * Rtot - rv * rv)
            inv_r = -1.0 / rv

            for iz in range(1, Zgrid - 1):
                for ku in range(num_u):
                    k = u[ku]

                    # -- diffusion --
                    pa = inv_r * (initc[ir, iz, k] - initc[ir - 1, iz, k]) / dr
                    pb = (initc[ir + 1, iz, k] - 2.0 * initc[ir, iz, k]
                          + initc[ir - 1, iz, k]) / dr2
                    # Non-uniform grid second derivative for axial diffusion
                    dx_c = dx[ir, iz, k]
                    dx_f = dx[ir, iz + 1, k] if iz + 1 < Zgrid else dx_c
                    dx_b = dx[ir, iz - 1, k] if iz - 1 >= 0 else dx_c
                    pc = (2.0 * (initc[ir, iz + 1, k] - initc[ir, iz, k])
                          / (dx_c * (dx_f + dx_c))
                          - 2.0 * (initc[ir, iz, k] - initc[ir, iz - 1, k])
                          / (dx_c * (dx_c + dx_b)))
                    diff = D[ir, iz, k] * (pa + pb + pc)

                    # -- convection --
                    adv = flow_fac * (initc[ir, iz, k] - initc[ir, iz - 1, k]) / dx_c

                    c[ir, iz, k] = dt * (diff - adv) + initc[ir, iz, k]

            # --- last column (iz = Zgrid-1) ---
            iz = Zgrid - 1
            for ku in range(num_u):
                k = u[ku]
                pa_e = inv_r * (initc[ir, iz, k] - initc[ir - 1, iz, k]) / dr
                pb_e = (initc[ir + 1, iz, k] - 2.0 * initc[ir, iz, k]
                        + initc[ir - 1, iz, k]) / dr2
                # End-column axial diffusion: zero-gradient outlet BC
                dx_val = dx[ir, iz, k]
                dx_b = dx[ir, iz - 1, k]
                pc_e = -2.0 * (initc[ir, iz, k] - initc[ir, iz - 1, k]) / (dx_val * (dx_val + dx_b))
                diff_e = D[ir, iz, k] * (pa_e + pb_e + pc_e)
                adv_e = flow_fac * (initc[ir, iz, k] - initc[ir, iz - 1, k]) / dx_val
                c[ir, iz, k] = dt * (diff_e - adv_e) + initc[ir, iz, k]

        # Row 0: term1[0]=term2[0]=0, so c[0] = 0 + initc[0] = unchanged
        # (already the case since we only write c[1:sp] above)

        # Symmetry: mirror upper half
        for ir in range(sp, Rgrid):
            mirror = Rgrid - 1 - ir
            for iz in range(Zgrid):
                for ku in range(num_u):
                    k = u[ku]
                    c[ir, iz, k] = c[mirror, iz, k]

    return c


def odesolve(timesteps, Zgrid, Rgrid, dt, D, Rtot, dr, dx, Qtot, c, u, rrc, modelparams, OHsource):
    """PDE solver for diffusion + convection in cylindrical tube.

    For model_mode='flowtube2', chemistry is handled externally by ode_solv_batch,
    so only diffusion and convection are computed here.
    For model_mode='flowtube1', chemistry term3 is computed in pure Python (fallback).
    """
    num = len(modelparams.comp_namelist)
    sp = int(Rgrid) // 2
    is_flowtube1 = (modelparams.model_mode == 'flowtube1')

    if not is_flowtube1:
        # flowtube2 / kinetic: pure diffusion + convection (no inline chemistry)
        # Ensure arrays are float64 contiguous for Numba
        c_f  = np.ascontiguousarray(c, dtype=np.float64)
        D_f  = np.ascontiguousarray(D, dtype=np.float64)
        dx_f = np.ascontiguousarray(dx, dtype=np.float64)
        u_a  = np.asarray(u, dtype=np.int64)

        c_f = _odesolve_inner(int(timesteps), int(Zgrid), int(Rgrid),
                              float(dt), D_f, float(Rtot), float(dr), dx_f,
                              float(Qtot), c_f, u_a, sp)
        return c_f
    else:
        # flowtube1 mode: fall back to original vectorized Python with inline chemistry
        initc = c  # alias, same as original
        _r1d = np.abs((np.arange(int(Rgrid), dtype=np.float64) - int(Rgrid) / 2 + 0.5) * float(dr))
        r = np.broadcast_to(_r1d[:, np.newaxis, np.newaxis], (int(Rgrid), int(Zgrid), num)).copy()

        term1 = np.zeros([int(Rgrid), int(Zgrid), num])
        term2 = np.zeros([int(Rgrid), int(Zgrid), num])
        term3 = np.zeros([int(Rgrid), int(Zgrid), num])

        dr2 = float(dr) ** 2

        for m in range(timesteps):
            p_a = -1. / r[1:sp, 1:-1, u] * (initc[1:sp, 1:-1, u] - initc[0:sp - 1, 1:-1, u]) / float(dr)
            p_b = (initc[2:sp + 1, 1:-1, u] - 2. * initc[1:sp, 1:-1, u] + initc[0:sp - 1, 1:-1, u]) / dr2
            p_c = (2.0 * (initc[1:sp, 2:, u] - initc[1:sp, 1:-1, u])
                   / (dx[1:sp, 1:-1, u] * (dx[1:sp, 2:, u] + dx[1:sp, 1:-1, u]))
                   - 2.0 * (initc[1:sp, 1:-1, u] - initc[1:sp, 0:-2, u])
                   / (dx[1:sp, 1:-1, u] * (dx[1:sp, 1:-1, u] + dx[1:sp, 0:-2, u])))
            term1[1:sp, 1:-1, u] = D[1:sp, 1:-1, u] * (p_a + p_b + p_c)

            conv_fac = (2. * Qtot) / (np.pi * Rtot ** 4) * (Rtot ** 2 - r[1:sp, 1:-1, u] ** 2)
            term2[1:sp, 1:-1, u] = conv_fac * (initc[1:sp, 1:-1, u] - initc[1:sp, 0:-2, u]) / dx[1:sp, 1:-1, u]

            p_a_end = -1. / r[1:sp, -1, u] * (initc[1:sp, -1, u] - initc[0:sp - 1, -1, u]) / float(dr)
            p_b_end = (initc[2:sp + 1, -1, u] - 2. * initc[1:sp, -1, u] + initc[0:sp - 1, -1, u]) / dr2
            p_c_end = -2.0 * (initc[1:sp, -1, u] - initc[1:sp, -2, u]) / (dx[1:sp, -1, u] * (dx[1:sp, -1, u] + dx[1:sp, -2, u]))
            term1[1:sp, -1, u] = D[1:sp, -1, u] * (p_a_end + p_b_end + p_c_end)

            conv_fac_end = (2. * Qtot) / (np.pi * Rtot ** 4) * (Rtot ** 2 - r[1:sp, -1, u] ** 2)
            term2[1:sp, -1, u] = conv_fac_end * (initc[1:sp, -1, u] - initc[1:sp, -2, u]) / dx[1:sp, -1, u]

            # Chemistry term for flowtube1 mode
            term3 = np.zeros([int(Rgrid), int(Zgrid), num])
            const_comp_set = set(modelparams.const_comp)
            for comp_na in modelparams.comp_namelist:
                if comp_na not in const_comp_set:
                    key_name = str(comp_na) + '_comp_indx'
                    compi = modelparams.dydt_vst[key_name]
                    key_name = str(comp_na) + '_res'
                    dydt_rec = modelparams.dydt_vst[key_name]
                    key_name = str(comp_na) + '_reac_sign'
                    reac_sign = modelparams.dydt_vst[key_name]
                    reac_count = 0
                    for i in dydt_rec[0, :]:
                        i = int(i)
                        if modelparams.Init_set == 'on':
                            gprate = initc[1:sp, 1:, modelparams.rindx_g[i, 0:modelparams.nreac_g[i]]] ** modelparams.rstoi_g[i, 0:modelparams.nreac_g[i]]
                        else:
                            gprate = initc[1:sp, 0:, modelparams.rindx_g[i, 0:modelparams.nreac_g[i]]] ** modelparams.rstoi_g[i, 0:modelparams.nreac_g[i]]
                        if len(modelparams.rstoi_g[i, 0:modelparams.nreac_g[i]]) > 1:
                            gprate1 = gprate[:, :, 0] * gprate[:, :, -1] * rrc[i]
                        else:
                            gprate1 = gprate[:, :, 0] * rrc[i]
                        if modelparams.Init_set == 'on':
                            term3[1:sp, 1:, compi] += reac_sign[reac_count] * gprate1
                        else:
                            term3[1:sp, 0:, compi] += reac_sign[reac_count] * gprate1
                        reac_count += 1

            c[0:sp, :, u] = dt * (term1[0:sp, :, u] - term2[0:sp, :, u] + term3[0:sp, :, u]) + initc[0:sp, :, u]
            c[sp:, :, u] = np.flipud(c[0:sp, :, u])
            initc = c

        return c