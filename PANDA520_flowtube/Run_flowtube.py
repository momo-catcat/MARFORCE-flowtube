## This file is used to run the flowtube model
import os
import numpy as np
import pandas as pd
import csv
from cmd_calib5 import cmd_calib5
from Calcu_by_flow import const_comp_conc_cal, const_comp_conc_cal_H2O, const_comp_conc_cal_OH
from diffusion_const_added import add_diff_const as add_diff_const


# add unit after values
class UnitFloat(float):

    def __new__(self, value, unit=None):
        return float.__new__(self, value)

    def __init__(self, value, unit=None):
        self.unit = unit


def Run_flowtube(para,folder_flowtube,export_file_folder,inputs_setup):

    for i in range(len(para['date'])):
        R1, L1, R2, L2, flag_tube, file, H2O_1, H2O_2, H2Oconc_1, H2Oconc_2, Q1, Q2, Itx, Qx, SO2flow, O2flow, N2flow, T_cel = inputs_setup(para['date'][i])

        T_cel=np.mean(T_cel)
        T = T_cel+273.15 # Temperature in K
        p = para['P'][i] # Pressure in Pa

        outflowLocation = para['outflowLocation'][i]  # outflow tube located before or after injecting air, water, and so2

        fullOrSimpleModel = para['fullOrSimpleModel'][i]  # simple: Gormley&Kennedy approximation, full: flow model (much slower)
        # in this study, we have the outflow before injecting air, water and SO2
        sampflow = para['sampleflow'][i]  # lpm
        O2ratio = para['O2ratio'][i]  # O2inAir = 0.209
        SO2ratio = para['SO2ratio'][i]

        const_comp_pre = ['SO2', 'O2']  # species have constant concentration and are calculated from flows
        const_comp_pre_know = ['H2O']  # species have known constant concentration but already known
        const_comp = const_comp_pre + const_comp_pre_know  # species have constant concentration
        # get all the concentrations
        O2conc, SO2conc = const_comp_conc_cal(O2flow, SO2flow,outflowLocation, sampflow, H2O_1,  N2flow, O2ratio, SO2ratio,
                                    Q1, Q2, T_cel, T, p, flag_tube)


        H2Oconc = const_comp_conc_cal_H2O(O2flow, SO2flow, outflowLocation, sampflow, H2O_1, H2O_2, N2flow, O2ratio,
                                        Q1, Q2, T_cel, T, p, flag_tube)


        # % store all the const species to const_comp_conc follow the order of const_comp
        const_comp_conc = np.transpose([SO2conc, O2conc, H2Oconc])

        OHconc, const_comp_free, const_comp_conc_free = const_comp_conc_cal_OH(H2Oconc, O2conc, Q1, flag_tube, Itx, Qx)

        # store initial concentration
        if para['model_mode'] == 'kinetic': # This mode runs the code in kinetic mode in which chemistry does not exist
            # define the H2SO4 concentration as 1e8 for convenience. !!!! This needs improvement
            OHconc = OHconc * 0 + 1e8
            Init_comp = ['OH', 'HO2', 'H2SO4']
            Init_comp_conc = np.transpose([OHconc, OHconc, OHconc])
        else:
            Init_comp = ['OH', 'HO2']  # species have initial concentration
            Init_comp_conc = np.transpose([OHconc, OHconc])


        print(Init_comp_conc)
        print(Init_comp)
        # add some diffusion constants add more in diffusion_const_added.py file if you want
        # set diffusivity according Kurten et al. (10.1021/jp993622j)
        Diff_setname = ['OH', 'H2O', 'HO2', 'SO3', 'H2SO4'] # the H2SO4 diffusivity is not final
        dH2O = add_diff_const(p, T)
        Diff_set = [0.215, dH2O, 0.141, 0.126, 0.088]

        # grid for length and radius directions
        Zgrid = np.array(para['Zgrid_num'][i]).astype(int)  # number of grid points in tube length direction
        Rgrid = np.array(para['Rgrid_num'][i]).astype(int)  # number of grid points in tube radius direction

        # chemistry part
        sch_name = folder_flowtube + 'input_mechanism/SO2_SA.txt'
        chm_sch_mrk = ['{', 'RO2', '+', '', '', ';', '+', ';', '$', '{', ':', ';', '}']
        key_spe_for_plot = 'H2SO4'
        plot_spec = ['OH', 'HSO3', 'HO2', 'SO3', 'H2SO4']  # plot species

        params = {'T': UnitFloat(T, "K"),  # temperature
                'p': UnitFloat(p, "Pa"),  # pressure pa
                'R1': UnitFloat(R1, "cm"),  # diameters for first tube
                'R2': UnitFloat(R2, "cm"),  # diameters for second tube
                'L1': UnitFloat(L1, "cm"),  # length for first tube
                'L2': UnitFloat(L2, "cm"),  # length for first tube
                'dt': para['dt'][i],  # dt * timesteps * numLoop is time elapsed in the final solution
                'Diff_setname': Diff_setname,  # diffusion for the species that you want to have
                'Diff_set': Diff_set,
                'fullOrSimpleModel': fullOrSimpleModel,  # Gormley&Kennedy approximation, full: flow model (much slower)
                'sch_name': sch_name,  # file for the MCM file
                'chm_sch_mrk': chm_sch_mrk,  # markers to isolate sections of chemical scheme based on MCM KPP format
                'const_comp': const_comp,  # constant species
                'Init_comp': Init_comp,  # species have initial concentration
                'Zgrid': Zgrid,  # number of grid points in tube length direction
                'Rgrid': Rgrid,  # number of grid points in tube radius direction
                # 'formula': formula, # the formula for the plots
                'key_spe_for_plot': key_spe_for_plot,  # key species for plot and criterion to stop the loop
                'plot_spec': plot_spec,  # plot species
                'flag_tube': flag_tube,
                'const_comp_free': const_comp_free,
                'const_comp_conc_free': const_comp_conc_free,
                'model_mode': para['model_mode']
                }
        # %

        #print('H2O',H2Oconc)
        #print('OH', OHconc)
        #print('const_comp_conc', const_comp_conc[:,1,:])
        #print('Init_comp_conc', Init_comp_conc[1])
        #print('params',params)
        #print('Q1',Q1[1])
        #print('Q2',Q2[1])

        #### SA






        # %% computation begins
        meanconc = []

        c = []

        for j in range(len(H2O_1)):
            if OHconc[j] > 0:
                meanConc1, c1 = cmd_calib5(const_comp_conc[:, j, :], params, Init_comp_conc[j], Q1[j], Q2[j])
                meanconc.append(meanConc1)
                print(meanconc)
                c.append(c1)

        meanconc_s = pd.DataFrame(meanconc,columns=plot_spec)


        # % save the modelled SA, HO2
        meanconc_s.to_csv(export_file_folder + 'SA_cali_' + str(para['date'][i]) + '.csv')

        with open(export_file_folder + 'SA_cali_' + str(para['date'][i]) + '.txt', 'w') as f:
            # using csv.writer method from CSV package
            write = csv.writer(f)

            write.writerows(c)
