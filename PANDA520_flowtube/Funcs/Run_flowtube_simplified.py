## This file is used to run the flowtube model

import pandas as pd
import csv
from Funcs.cmd_calib5 import cmd_calib5


#%%
def Run_flowtube(paras, export_file_folder, const_comp_conc, Init_comp_conc, num_stage):
    #%%
    meanconc = []
    c = []
    if isinstance(num_stage, int):
        for j in range(num_stage):
            meanConc1, c1 = cmd_calib5(const_comp_conc[:, j, :], paras, Init_comp_conc[j], paras['Q1'][j],
                                       paras['Q2'][j])
            meanconc.append(meanConc1)
            c.append(c1)
    else: # SA calibraiton case
        for j in range(len(num_stage)):
            if num_stage[j] > 0:
                print(j)
                meanConc1, c1 = cmd_calib5(const_comp_conc[:, j, :], paras, Init_comp_conc[j], paras['Q1'][j],paras['Q2'][j])
                #plt.close()
                meanconc.append(meanConc1)
                c.append(c1)


    meanconc_s = pd.DataFrame(meanconc, columns=paras['plot_spec'])
    meanconc_s.to_csv(export_file_folder + str(paras['file_name'][:-4]) + '_output.csv')

    with open(export_file_folder + str(paras['file_name'][:-4]) + '_output.txt', 'w') as f:
        # using csv.writer method from CSV package
        write = csv.writer(f)
        write.writerows(c)
