'''function to initiate concentrations of components'''
# #########################################################################################
# Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk
# # All Rights Reserved.
# This file is part of PyCHAM
# # PyCHAM is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
# # PyCHAM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more details.
# # You should have received a copy of the GNU General Public License along with PyCHAM.
# If not, see <http://www.gnu.org/licenses/>.
# #########################################################################################
# based on inputs, initial concentrations and their holding arrays 
# are set

import numpy as np
import scipy.constants as si
from Funcs import write_dydt_rec


def init_conc(num_comp, Comp0, init_conc, TEMP, PInit,
              dydt_trak, rindx, pindx, num_eqn, nreac, nprod,
              comp_namelist, RO2_indx, HOMRO2_indx, rstoi, pstoi):
    # inputs:------------------------------------------------------

    # num_comp - number of unique components
    # Comp0 - chemical scheme names of components present at start of experiment
    # init_conc - initial concentrations of components (ppb)
    # TEMP - temperature in chamber at start of experiment (K)
    # PInit - initial pressure (Pa)
    # dydt_trak - chemical scheme name of components for which user wants the tendency to  change tracked
    # end_sim_time - total simulation time (s)
    # save_step - recording frequency (s)
    # rindx - indices of reactants per equation
    # pindx - indices of products per equation
    # num_eqn - number of equations
    # comp_namelist - list of names of components as presented in the chemical scheme file
    # nuc_comp - name of nucleating component (input by user, or defaults to 'core')
    # RO2_indx - RO2 list indices and chemical scheme indices of non-HOM-RO2 molecules
    # HOMRO2_indx - chemical scheme indices of HOM-RO2 molecules
    # rstoi - stoichiometry of reactants per equation
    # pstoi - stoichiometry of products per equation
    # -----------------------------------------------------------

    # start by assuming no error
    erf = 0
    err_mess = ''

    # 	if testf==1: # testing mode
    # 		# return dummies
    # 		return(0,0,0,0,0,0,0,0)

    NA = si.Avogadro  # Avogadro's number (molecules/mol)
    # empty array for storing species' concentrations, must be an array
    y = np.zeros((num_comp))
    y_mw = np.zeros((num_comp, 1))  # species' molecular weight (g/mol)
    # empty array for storing index of interesting gas-phase components
    y_indx_plot = []

    # convert concentrations
    # total number of molecules in 1 cc air using ideal gas law.  R has units cc.Pa/K.mol
    ntot = PInit * (NA / ((si.R * 1.e6) * TEMP))
    # one billionth of number of # molecules in chamber unit volume
    Cfactor = ntot * 1.e-9  # ppb to # molecules/cm3 conversion factor

    # prepare dictionary for tracking tendency to change of user-specified components
    dydt_vst = {}

    # insert initial concentrations where appropriate
    for i in range(len(Comp0)):
        # index of where initial components occur in list of components
        try:  # in case components already listed via interpretation of the chemical scheme
            y_indx = comp_namelist.index(Comp0[i])

        # if component not already listed via interpretation of the chemical scheme
        # then send error message
        except:
            erf = 1
            err_mess = str('Error: component called ' + str(Comp0[
                                                                i]) + ', which has an initial concentration specified in the model variables input file has not been found in the chemical scheme.  Please check the scheme and associated chemical scheme markers, which are stated in the model variables input file.')
            return (0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0,
                    0, 0, 0, 0, erf, err_mess, 0, 0, 0)

        y[y_indx] = init_conc[i]  # *Cfactor # convert from ppb to # molecules/cm3 (air)

        # remember index for plotting gas-phase concentrations later
        y_indx_plot.append(y_indx)

    # --------------------------------------
    # get index of user-specified components for tracking their change tendencies (dydt) due to modeled
    # mechanisms
    if (len(dydt_trak) > 0):

        dydt_traki = []  # empty list for indices of these components

        for i in range(len(dydt_trak)):

            reac_index = []  # indices of reactions involving this component
            # value for whether this component is reactant or product in a reaction and stoichiometry
            reac_sign = []

            if dydt_trak[i] != 'RO2' and dydt_trak[i] != 'HOMRO2':

                # index of components in component list
                try:
                    y_indx = comp_namelist.index(dydt_trak[i])
                # if component not already listed via interpretation of the chemical scheme
                # then send error message
                except:
                    erf = 1
                    err_mess = str('Error: component called ' + str(dydt_trak[
                                                                        i]) + ', which is specified to be tracked in the model variables input file has not been found in the chemical scheme.  Please check the scheme and associated chemical scheme markers, which are stated in the model variables input file.')
                    return (0, 0, 0, 0, 0, 0, 0, 0,
                            0, 0, 0,
                            0, 0, 0, 0, erf, err_mess, 0, 0, 0)
                # remember index for plotting gas-phase concentrations later
                dydt_traki.append([int(y_indx)])

                # search through reactions to see where this component is reactant or product
                for ri in range(num_eqn):
                    if sum(rindx[ri, 0:nreac[ri]] == y_indx) > 0:
                        reac_index.append(int(ri))  # append reaction index
                        reac_place = np.where(rindx[ri, 0:nreac[ri]] == y_indx)[0]
                        reac_sign.append(-1 * rstoi[int(ri), reac_place])
                    if sum(pindx[ri, 0:nprod[ri]] == y_indx) > 0:
                        reac_index.append(int(ri))  # append reaction index
                        reac_place = np.where(pindx[ri, 0:nprod[ri]] == y_indx)[0]
                        reac_sign.append(1 * pstoi[int(ri), reac_place])

            if dydt_trak[i] == 'RO2':  # for non-HOM-RO2 (organic peroxy radicals)

                # remember index for plotting gas-phase concentrations later
                dydt_traki.append(list(RO2_indx[:, 1]))

                # loop through non-HOM-RO2 components to get their reaction indices
                for y_indx in RO2_indx[:, 1]:

                    # search through reactions to see where this component is reactant or product
                    for ri in range(num_eqn):
                        if sum(rindx[ri, 0:nreac[ri]] == y_indx) > 0:
                            reac_index.append(int(ri))  # append reaction index
                            reac_place = np.where(rindx[ri, 0:nreac[ri]] == y_indx)[0]
                            reac_sign.append(-1 * rstoi[int(ri), reac_place])

                        if sum(pindx[ri, 0:nprod[ri]] == y_indx) > 0:
                            reac_index.append(int(ri))  # append reaction index
                            reac_place = np.where(pindx[ri, 0:nprod[ri]] == y_indx)[0]
                            reac_sign.append(1 * pstoi[int(ri), reac_place])

                    y_indx = RO2_indx[:, 1]  # ready for storing below

            if dydt_trak[i] == 'HOMRO2':  # for HOM-RO2 (highly oxidised molecule organic peroxy radicals)

                # remember index for plotting gas-phase concentrations later
                dydt_traki.append(list(np.squeeze(HOMRO2_indx)))

                # loop through non-HOM-RO2 components to get their reaction indices
                for y_indx in HOMRO2_indx[:]:

                    # search through reactions to see where this component is reactant or product
                    for ri in range(num_eqn):
                        if sum(rindx[ri, 0:nreac[ri]] == y_indx) > 0:
                            reac_index.append(int(ri))  # append reaction index
                            reac_place = np.where(rindx[ri, 0:nreac[ri]] == y_indx)[0]
                            reac_sign.append(-1 * rstoi[int(ri), reac_place])
                        if sum(pindx[ri, 0:nprod[ri]] == y_indx) > 0:
                            reac_index.append(int(ri))  # append reaction index
                            reac_place = np.where(pindx[ri, 0:nprod[ri]] == y_indx)[0]
                            reac_sign.append(1 * pstoi[int(ri), reac_place])
                y_indx = HOMRO2_indx[:]  # ready for storing below

            # save reaction indices in dictionary value for this component,
            # when creating empty rec_array, add two columns onto the end for particle- and
            # wall-partitioning, respectively
            rec_array = np.zeros((2, len(reac_index)))
            rec_array[0, :] = reac_index

            comp_indx_str = str(dydt_trak[i] + '_comp_indx')
            res_string = str(dydt_trak[i] + '_res')
            reac_string = str(dydt_trak[i] + '_reac_sign')

            dydt_vst[comp_indx_str] = y_indx  # dictionary entry to hold index of tracked component
            dydt_vst[res_string] = rec_array  # dictionary entry to hold reaction indices and results
            dydt_vst[reac_string] = reac_sign  # dictionary entry to hold results

        # dictionary entry to hold component names of components to track
        dydt_vst['comp_names'] = dydt_trak

        # call on write_dydt_rec to generate the module that will process
        # the tendency to change during the simulation
        write_dydt_rec.write_dydt_rec()

    return y, Cfactor * 1e9, dydt_vst
