def odesolve(timesteps, Zgrid, Rgrid, dt, D, Rtot, dr, dx, Qtot, c, comp_namelist, dydt_vst, rindx, nreac, rstoi,
             rate_values, const_comp, u, model_mode):
    # %import packages
    import numpy as np
    # %%
    # D = Diff_vals
    initc = c
    num = len(comp_namelist)
    r = np.zeros([int(Rgrid), int(Zgrid), num])
    r = np.abs(np.array([r[i, :, :] + i - r.shape[0] / 2 + 0.5 for i in range(r.shape[0])])) * dr

    # do calculations in matrix
    D_temp = np.tile(D, (Rgrid * Zgrid))
    D = np.transpose(np.reshape(D_temp, (Rgrid, Zgrid, num)), (0, 1, 2))

    # define production and loss terms
    term1 = np.zeros([int(Rgrid), int(Zgrid), num])
    term2 = np.zeros([int(Rgrid), int(Zgrid), num])
    term3 = np.zeros([int(Rgrid), int(Zgrid), num])

    # timesteps = 173
    # %%
    for m in range(timesteps):
        # Diffusion term
        # The equation is based on Fick's first law: https://en.wikipedia.org/wiki/Fick%27s_laws_of_diffusion
        # calculate central grids

        p_a = - 1. / r[1:Rgrid // 2, 1:-1, u] * (initc[1:Rgrid // 2, 1:-1, u] - initc[0:Rgrid // 2 - 1, 1:-1, u]) / dr[
                                                                                                                    1:Rgrid // 2,
                                                                                                                    1:-1,
                                                                                                                    u]

        p_b = (initc[2:Rgrid // 2 + 1, 1:-1, u] - 2. * initc[1:Rgrid // 2, 1:-1, u] + initc[0:Rgrid // 2 - 1, 1:-1,
                                                                                      u]) / (
                      dr[1:Rgrid // 2, 1:-1, u] ** 2)

        p_c = (initc[1:Rgrid // 2, 2:, u] - 2. * initc[1:Rgrid // 2, 1:-1, u] + initc[1:Rgrid // 2, 0:-2, u]) / (
                dx * dx)

        term1[1:Rgrid // 2, 1:-1, u] = D[1:Rgrid // 2, 1:-1, u] * (p_a + p_b + p_c)  # diffusion of gas molecules

        # convection; carried by main flow
        # Refs: 1. https://en.wikipedia.org/wiki/Advection
        #       2. Gormley & Kennedy, 1948, Diffusion from a stream flowing through a cylindrical tube
        term2[1:Rgrid // 2, 1:-1, u] = (2. * Qtot) / (np.pi * Rtot[1:Rgrid // 2, 1:-1, u] ** 4) * \
                                       (Rtot[1:Rgrid // 2, 1:-1, u] ** 2 - r[1:Rgrid // 2, 1:-1, u] ** 2) * \
                                       (initc[1:Rgrid // 2, 1:-1, u] - initc[1:Rgrid // 2, 0:-2,
                                                                       u]) / dx  # carried by main flow

        # calculate the last column (measured by the instrument)

        p_a_end = - 1. / r[1:Rgrid // 2, -1, u] * (
                initc[1:Rgrid // 2, -1, u] - initc[0:Rgrid // 2 - 1, -1, u]) / dr[1:Rgrid // 2, -1, u]

        p_b_end = (initc[2:Rgrid // 2 + 1, -1, u] - 2. * initc[1:Rgrid // 2, -1, u] + initc[0:Rgrid // 2 - 1, -1,
                                                                                      u]) / (
                          dr[1:Rgrid // 2, -1, u] ** 2)

        p_c_end = (initc[1:Rgrid // 2, -1, u] - 2. * initc[1:Rgrid // 2, -2, u] + initc[1:Rgrid // 2, -2, u]) / (
                dx * dx)

        term1[1:Rgrid // 2, -1, u] = D[1:Rgrid // 2, -1, u] * (p_a_end + p_b_end + p_c_end)

        term2[1:Rgrid // 2, -1, u] = (2. * Qtot) / (np.pi * Rtot[1:Rgrid // 2, -1, u] ** 4) * \
                                     (Rtot[1:Rgrid // 2, -1, u] ** 2 - r[1:Rgrid // 2, -1, u] ** 2) * \
                                     (initc[1:Rgrid // 2, -1, u] - initc[1:Rgrid // 2, -2,
                                                                   u]) / dx  # carried by main flow



        if model_mode == 'kinetic':
            c[0:Rgrid // 2, :, u] = dt * (
                    term1[0:Rgrid // 2, :, u] - term2[0:Rgrid // 2, :, u]) + initc[0:Rgrid // 2,:, u]

        else:
            term3 = np.zeros([int(Rgrid), int(Zgrid), num])

            for comp_na in comp_namelist:  # get name of this component
                #if comp_na not in const_comp:
                key_name = str(str(comp_na) + '_comp_indx')  # get index of this component
                compi = dydt_vst[key_name]
                key_name = str(str(comp_na) + '_res')
                dydt_rec = dydt_vst[key_name]
                key_name = str(str(comp_na) + '_reac_sign')
                reac_sign = dydt_vst[key_name]
                # dydt_for = np.zeros(comp_num)
                reac_count = 0
                for i in dydt_rec[0, :]:
                    i = int(
                            )  # ensure reaction index is integer - this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float)
                    gprate = initc[1:Rgrid // 2, 1:, rindx[i, 0:nreac[i]]] ** rstoi[i, 0:nreac[i]]
                    if (len(rstoi[i, 0:nreac[i]]) > 1):

                        # print(rate_values)
                        # print(gprate[:,:,0], gprate[:, :, -1], rate_values[i])
                        gprate1 = gprate[:, :, 0] * gprate[:, :, -1] * rate_values[i]
                        # gprate2 = gprate[:,:,0] * gprate[:, :,-1] * rate_values[i]
                        # gprate3 = gprate1-gprate2
                    else:
                        gprate1 = gprate[:, :, 0] * rate_values[i]

                    term3[1:Rgrid // 2, 1:, compi] += reac_sign[reac_count] * (gprate1)
                    reac_count += 1

            c[0:Rgrid // 2, :, u] = dt * (
                    term1[0:Rgrid // 2, :, u] - term2[0:Rgrid // 2, :, u] + term3[0:Rgrid // 2, :, u]) + initc[
                                                                                                         0:Rgrid // 2,
                                                                                                         :, u]
        #c[0:Rgrid // 2, :, u] = dt * (
        #            term1[0:Rgrid // 2, :, u] - term2[0:Rgrid // 2, :, u] ) + initc[
        #                                                                                                 0:Rgrid // 2,
        #                                                                                                 :, u]

        # c[0:Rgrid // 2, :, u] = dt * (term1[0:Rgrid // 2, :, u] - term2[0:Rgrid // 2, :, u])  + initc[0:Rgrid // 2,:, u]

        c[Rgrid // 2:, :, u] = np.flipud(c[0:Rgrid // 2, :, u])

        initc = c
    # %%
    return (c)