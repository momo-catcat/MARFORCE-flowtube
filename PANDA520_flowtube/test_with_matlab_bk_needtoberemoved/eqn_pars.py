'''parses the input files to automatically create the solver file'''
# #########################################################################################
# Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk
# # All Rights Reserved.
# This file is part of PyCHAM, here, we have small modification
# # PyCHAM is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
# # PyCHAM is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more details.
# # You should have received a copy of the GNU General Public License along with PyCHAM.
# If not, see <http://www.gnu.org/licenses/>.
# #########################################################################################
# input files are interpreted and used to create the necessary
# arrays and python files to solve problem

import numpy as np
from Funcs import eqn_interr, RO2_indices, sch_interr, write_dydt_rec


# define function to extract the chemical mechanism
def extr_mech(sch_name, chem_sch_mrk,
              con_infl_nam, const_comp):
    #          drh_str, erh_str):
    # inputs: ----------------------------------------------------
    # sch_name - file name of chemical scheme
    # chem_sch_mrk - markers to identify different sections of
    # 	the chemical scheme
    # con_infl_nam - chemical scheme names of components with
    # 		constant influx
    # int_tol - integration tolerances
    # wall_on - marker for whether to include wall partitioning
    # num_sb - number of size bins (including any wall)
    # const_comp - chemical scheme name of components with
    #	constant concentration
    # ------------------------------------------------------------

    # starting error flag and message (assumes no errors)
    erf = 0
    err_mess = ''

    f_open_eqn = open(sch_name, mode='r')  # open the chemical scheme file
    # read the file and store everything into a list
    total_list_eqn = f_open_eqn.readlines()
    f_open_eqn.close()  # close file

    # interrogate scheme to list equations
    [eqn_list, eqn_num, rrc, rrc_name,
     RO2_names, eqn_list_on] = sch_interr.sch_interr(chem_sch_mrk, sch_name)

    # interrogate xml to list all component names and SMILES
    # 	[comp_smil, comp_name] = xml_interr.xml_interr(xml_name)

    # get equation information for chemical reactions
    [rindx_g, rstoi_g, pindx_g, pstoi_g, reac_coef_g,
     nreac_g, nprod_g, comp_namelist, comp_num] = eqn_interr.eqn_interr(eqn_num,
                                                                        eqn_list, chem_sch_mrk)

  # get index of components with constant influx/concentration -----------
    # empty array for storing index of components with constant influx
    con_infl_indx = np.zeros((len(con_infl_nam)))
    con_C_indx = np.zeros((len(const_comp))).astype('int')
    for i in range(len(con_infl_nam)):
        # water not included explicitly in chemical schemes but accounted for later in init_conc
        if (con_infl_nam[i] == 'H2O'):
            con_infl_indx[i] = comp_num
            continue
        try:
            # index of where components with constant influx occur in list of components
            con_infl_indx[i] = comp_namelist.index(con_infl_nam[i])
        except:
            erf = 1  # raise error
            err_mess = str('Error: constant influx component with name ' + str(con_infl_nam[
                                                                                   i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')

    for i in range(len(const_comp)):
        try:
            # index of where constant concentration components occur in list
            # of components
            con_C_indx[i] = comp_namelist.index(const_comp[i])
        except:
            erf = 1  # raise error
            err_mess = str('Error: constant concentration component with name ' + str(const_comp[
                                                                                          i]) + ' has not been identified in the chemical scheme, please check it is present and the chemical scheme markers are correct')

    # ---------------------------------------------------------------------

    # call function to generate ordinary differential equation (ODE)
    # solver module, add two to comp_num to account for water and core component
    # 	write_ode_solv.ode_gen(con_infl_indx, int_tol, rowvals, wall_on, comp_num+2,
    # 			(num_sb-wall_on), 0, eqn_num, dil_fac, sav_nam, pcont)

    # call function to generate reaction rate calculation module
    # 	write_rate_file.write_rate_file(reac_coef_g, reac_coef_aq, rrc, rrc_name, 0)

    # call function to generate module that tracks change tendencies
    # of certain components
    write_dydt_rec.write_dydt_rec()

    # write the module for estimating deliquescence and efflorescence
    # relative humidities as a function of temperature
    # write_hyst_eq.write_hyst_eq(drh_str, erh_str)

    # get index of components in the peroxy radical list
    RO2_indx = RO2_indices.RO2_indices(comp_namelist, RO2_names)

    # get index of HOM-RO2 radicals
    HOMRO2_indx = RO2_indices.HOMRO2_indices(comp_namelist)

    # get number of photolysis equations
    # Jlen = photo_num.photo_num(photo_path)

    return (RO2_indx, HOMRO2_indx, con_C_indx)
