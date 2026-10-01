## This file is used to run the flowtube model

import pandas as pd
import csv
from Funcs.cmd_calib5 import cmd_calib5
import numpy as np

#%%
def Run_flowtube(modelparams,num_stage):
    #%%
    meanconc = []
    c = []
    c_prev = None  # pass converged state from stage N to stage N+1

    if 'box' in modelparams.model_mode:
        # Box model: single stage, use first row of input data
        modelparams.number_stage = 0
        const_conc = modelparams.const_comp_conc[0] if modelparams.const_comp_conc.ndim >= 2 else modelparams.const_comp_conc
        if hasattr(modelparams, 'Init_comp_conc'):
            init_conc = modelparams.Init_comp_conc[0] if modelparams.Init_comp_conc.ndim >= 2 else modelparams.Init_comp_conc
        else:
            init_conc = np.array([], dtype=np.float32)
        q1 = modelparams.Q1[0] if hasattr(modelparams.Q1, '__len__') else modelparams.Q1
        q2 = modelparams.Q2[0] if hasattr(modelparams.Q2, '__len__') else modelparams.Q2
        # Reshape const_conc to (1, num_species) so cmd_calib5 can index [0, :]
        const_conc = const_conc.reshape(1, -1)
        meanConc1, c1 = cmd_calib5(const_conc, modelparams, init_conc, q1, q2, c_prev=None)
        meanconc.append(meanConc1)
        if isinstance(c1, list):
            c.append(c1[-1])
        else:
            c.append(c1)

    elif isinstance(num_stage, int):
        for j in range(num_stage):
            modelparams.number_stage = j
            meanConc1, c1 = cmd_calib5(modelparams.const_comp_conc[:, j, :], modelparams, modelparams.Init_comp_conc[j], modelparams.Q1[j],
                                       modelparams.Q2[j], c_prev=c_prev)
            meanconc.append(meanConc1)
            # c1 is (c, c_2nd) tuple for twotubes, c array for onetube, or list for box
            if isinstance(c1, tuple):
                c_out = c1[1]  # c_2nd for twotubes output
            elif isinstance(c1, list):
                c_out = c1[-1]  # last converged state for box model
            else:
                c_out = c1
            c.append(c_out)
            c_prev = c1
    else: # SA calibraiton case
        for j in range(len(num_stage)):
            modelparams.number_stage = j
            if num_stage[j] > 0:
                print(f"Running stage {j}")
                meanConc1, c1 = cmd_calib5(modelparams.const_comp_conc[:, j, :], modelparams, modelparams.Init_comp_conc[j], modelparams.Q1[j],modelparams.Q2[j], c_prev=c_prev)
                #plt.close()
                meanconc.append(meanConc1)
                # c1 is (c, c_2nd) tuple for twotubes, c array for onetube, or list for box
                if isinstance(c1, tuple):
                    c_out = c1[1]
                elif isinstance(c1, list):
                    c_out = c1[-1]
                else:
                    c_out = c1
                output_file = f"{modelparams.export_file_folder}R{modelparams.Rgrid}L{modelparams.Zgrid}inter_acc{modelparams.actual_acc:.0e}dt{modelparams.dt:.0e}timstep{modelparams.timesteps:.0e}itx{modelparams.Itx:.1e}{modelparams.model_mode}_Stage_{j}_endtube_profile_output.txt"
                with open(output_file, 'w', newline='', encoding='utf-8') as f:
                    write = csv.writer(f)
                    write.writerows(c_out)  # Write current dataset to a new file
                print(f"Saved stage {j + 1} to {output_file}")
                c_prev = c1

    meanconc_s = pd.DataFrame(meanconc, columns=modelparams.comp_namelist)
    meanconc_s.to_csv(f"{modelparams.export_file_folder}R{modelparams.Rgrid}L{modelparams.Zgrid}inter_acc{modelparams.actual_acc:.0e}dt{modelparams.dt:.0e}timstep{modelparams.timesteps:.0e}itx{modelparams.Itx:.1e}{modelparams.model_mode}_allstage_final_results.csv")

    exclude_keys = {"dydt_vst", "const_comp_gird"}  # Define keys to exclude

    with open(f"{modelparams.export_file_folder}R{modelparams.Rgrid}L{modelparams.Zgrid}inter_acc{modelparams.actual_acc:.0e}dt{modelparams.dt:.0e}timstep{modelparams.timesteps:.0e}itx{modelparams.Itx:.1e}{modelparams.model_mode}_params.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Parameter", "Value"])  # Add header

        for key, value in modelparams.__dict__.items():
            if key not in exclude_keys:
                writer.writerow([key, value])
