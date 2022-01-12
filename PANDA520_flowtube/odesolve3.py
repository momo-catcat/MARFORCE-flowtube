def odesolve(timesteps, Zgrid, Rgrid, dt, kSO2pOH, kOHpHO2, kOHpOH, kSO3p2H2O, kHSO3pO2,  D, Rtot, dr, dx, Qtot,c,comp_namelist,dydt_vst,rindx,nreac,rstoi,rate_values,SO2tot,H2Otot,O2tot):
    #%import packages
    import numpy as np
    import matplotlib.pyplot as plt
    initc = c
    num = 9
    r = np.zeros([int(Rgrid), int(Zgrid),  num])
    r = np.abs(np.array([r[i, :, :] + i - r.shape[0] / 2 + 0.5 for i in range(r.shape[0])])) * dr

    #do calculations in matrix
    D_temp = np.tile(D, (Rgrid * Zgrid))
    D = np.transpose(np.reshape(D_temp, (Rgrid, Zgrid, num)), (0, 1, 2))


    #define production and loss terms
    term1 = np.zeros([int(Rgrid), int(Zgrid), num])
    term2 = np.zeros([int(Rgrid), int(Zgrid), num])
    term3 = np.zeros([int(Rgrid), int(Zgrid), num])
    
    # SO2tot = c[:,:,comp_namelist.index('SO2')]
    # O2tot = c[:,:,comp_namelist.index('O2')]
    # H2Otot = c[:,:,comp_namelist.index('H2O')]

    const_comp = ['SO2','O2','H2O']
    u = [0,2,3,6,7,8]
    # u = [[0,1,2,3,4,5]

    for m in range(timesteps):
        # Diffusion term
        # The equation is based on Fick's first law: https://en.wikipedia.org/wiki/Fick%27s_laws_of_diffusion
        # calculate central grids
#%
        p_a = - 1. / r[1:Rgrid // 2, 1:-1, :] * (initc[1:Rgrid // 2, 1:-1, :] - initc[0:Rgrid // 2 - 1, 1:-1, :]) / dr[1:Rgrid // 2, 1:-1, :]

        p_b = (initc[2:Rgrid // 2 + 1, 1:-1, :] - 2. * initc[1:Rgrid // 2, 1:-1, :] + initc[0:Rgrid // 2 - 1, 1:-1, :]) / (dr[1:Rgrid // 2, 1:-1, :] ** 2)

        p_c = (initc[1:Rgrid // 2, 2:, :] - 2. * initc[1:Rgrid // 2, 1:-1, :] + initc[1:Rgrid // 2, 0:-2,: ]) / (dx * dx)

        term1[1:Rgrid // 2, 1:-1, :] = D[1:Rgrid // 2, 1:-1, :] * (p_a + p_b + p_c) #diffusion of gas molecules

        # convection; carried by main flow
        # Refs: 1. https://en.wikipedia.org/wiki/Advection
        #       2. Gormley & Kennedy, 1948, Diffusion from a stream flowing through a cylindrical tube
        term2[1:Rgrid // 2, 1:-1, :] = (2. * Qtot[1:Rgrid // 2, 1:-1,:]) / (np.pi * Rtot[1:Rgrid // 2, 1:-1,:] ** 4) * \
                                          (Rtot[1:Rgrid // 2, 1:-1,:] ** 2 - r[1:Rgrid // 2, 1:-1, :] ** 2) *\
                            (initc[1:Rgrid // 2, 1:-1, :] - initc[1:Rgrid // 2, 0:-2, :]) / dx # carried by main flow

        #calculate the last column (measured by the instrument)

        # p_a_end = - 1. / r[1:Rgrid // 2, -1, :] * (
        #             initc[1:Rgrid // 2, -1, :] - initc[0:Rgrid // 2 - 1, -1, :]) / dr[1:Rgrid // 2, -1, :]

        # p_b_end = (initc[2:Rgrid // 2 + 1, -1, :] - 2. * initc[1:Rgrid // 2, -1, :] + initc[0:Rgrid // 2 - 1, -1, :]) / (dr[1:Rgrid // 2, -1, :] ** 2)

        # p_c_end = (initc[1:Rgrid // 2, -1, :] - 2. * initc[1:Rgrid // 2, -2, :] + initc[1:Rgrid // 2, -2,: ]) / (dx * dx)

        # term1[1:Rgrid // 2, -1, :] = D[1:Rgrid // 2, -1, :] * (p_a_end + p_b_end + p_c_end)

        # term2[1:Rgrid // 2, -1, :] = (2. * Qtot[1:Rgrid // 2, -1,:]) / (np.pi * Rtot[1:Rgrid // 2, -1,:] ** 4) * \
        #                                 (Rtot[1:Rgrid // 2, -1,:] ** 2 - r[1:Rgrid // 2, -1, :] ** 2) *\
        #                     (initc[1:Rgrid // 2, -1, :] - initc[1:Rgrid // 2, -2, :]) / dx #carried by main flow

        # t = ['HSO_3', 'SO_3', 'HO_2', 'H_2SO_4', 'OH']
#         term3[1:-1, 1:, comp_namelist.index('HSO3')] = kSO2pOH *  SO2tot[1:-1, 1:] * initc[1:-1, 1:, comp_namelist.index('OH')] - kHSO3pO2 * initc[1:-1, 1:, comp_namelist.index('HSO3')] * O2tot[1:-1, 1:]

#         term3[1:-1, 1:, comp_namelist.index('SO3')] = kHSO3pO2 * O2tot[1:-1, 1:] * initc[1:-1, 1:, comp_namelist.index('HSO3')] - kSO3p2H2O * H2Otot[1:-1, 1:] * H2Otot[1:-1, 1:] * initc[1:-1, 1:, comp_namelist.index('SO3')]

#         term3[1:-1, 1:, comp_namelist.index('HO2')] = kHSO3pO2 * O2tot[1:-1, 1:] * initc[1:-1, 1:, comp_namelist.index('HSO3')] - kOHpHO2 * initc[1:-1, 1:, comp_namelist.index('HO2')] * initc[1:-1, 1:, comp_namelist.index('OH')]

#         term3[1:-1, 1:, comp_namelist.index('SA')] = kSO3p2H2O * H2Otot[1:-1, 1:] * H2Otot[1:-1, 1:] * initc[1:-1, 1:, comp_namelist.index('SO3')]

#         term3[1:-1, 1:,comp_namelist.index('OH')] = -kSO2pOH * SO2tot[1:-1, 1:] * initc[1:-1, 1:, comp_namelist.index('OH')] - kOHpHO2 * initc[1:-1, 1:, comp_namelist.index('HO2')] * initc[\
#             1:-1, 1:, comp_namelist.index('OH')] - 2. * kOHpOH * initc[1:-1, 1:, comp_namelist.index('OH')] * initc[1:-1, 1:, comp_namelist.index('OH')]
# # NEW VERSION
        term3 = np.zeros([int(Rgrid), int(Zgrid), num])


        for comp_na in comp_namelist: # get name of this component
            if comp_na not in const_comp: 
                key_name = str(str(comp_na)+ '_comp_indx')  # get index of this component
                compi = dydt_vst[key_name]
                key_name = str(str(comp_na)+'_res')
                dydt_rec = dydt_vst[key_name]
                key_name = str(str(comp_na) + '_reac_sign')
                reac_sign = dydt_vst[key_name]
                # dydt_for = np.zeros(comp_num)
                reac_count = 0
                for i in dydt_rec[0,:]:
                    i = int(i) # ensure reaction index is integer - this necessary because the dydt_rec array is float (the tendency to change records beneath its first row are float)
                    # print(i)
                    # print(rin1dx[i, 0:nreac[i]])
                    gprate = initc[1:Rgrid // 2, 1:,rindx[i, 0:nreac[i]]] ** rstoi[i, 0:nreac[i]]
                    if (len( rstoi[i, 0:nreac[i]]) > 1): 
                        gprate1= gprate[:,:,0] * gprate[:, :,-1] * rate_values[i]
                        
                    else:
                        gprate1= gprate[:,:,0]  * rate_values[i]
                    # print(gprate[5,2,:])

                    term3[1:Rgrid // 2, 1:, compi] += reac_sign[reac_count]*(gprate1)
                    # print(reac_sign[reac_count])
                    reac_count += 1
        # comp_na = 'OH'
        # key_name = str(str(comp_na)+ '_comp_indx')  # get index of this component
        # compi = dydt_vst[key_name]
        # oh1 = initc[:,:,compi]
        # ohso3 = initc[:,:,2]
        
        # hso31= term3[:,:,compi]
        # gp1 = c[1:-1, 1:,0]
        # gp11 = c[1:-1,1:,1]
        # gp2 = rstoi[i, 0:nreac[i]]
        # gp3 = (gp1**gp2[0])
        # gp33 = (gp11**gp2[1])
        # gp4= gp3[:,:]
        # gp5= gp3[:,:]
        # gp6 = gp3*gp33*rate_values[0]
        # term3[1:-1, 1:, compi] += reac_sign[reac_count]*((gp6))
# NEW VERSION        
        # c = dt * (term1 - term2 + term3) + initc
        c[0:Rgrid // 2,:,u] = dt * (term1[0:Rgrid // 2,:,u] - term2[0:Rgrid // 2,:,u] + term3[0:Rgrid // 2,:,u]) + initc[0:Rgrid // 2,:,u]
        #c[0:Rgrid // 2, :, u] = dt * (term1[0:Rgrid // 2, :, u] - term2[0:Rgrid // 2, :, u])  + initc[0:Rgrid // 2,:, u]

        c[Rgrid // 2:, :, u] = np.flipud(c[0:Rgrid // 2,:,u])
        # c[:,:,1] = SO2tot
        # c[:,:,5] = O2tot
        # c[:,:,4] = H2Otot
        initc = c
        oh1 = initc[:,:,0]
        ohso3 = initc[:,:,2]
        ohso3 = initc[:,:,7]
#% hso3 = c[:,:,comp_namelist.index('HSO3')]
    return(c)



