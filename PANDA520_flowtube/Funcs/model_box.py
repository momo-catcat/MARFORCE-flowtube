import numpy as np
from Funcs.meanconc_cal import meanconc_cal_sim as meanconc_cal_sim
from Funcs.odesolve3 import odesolve as odesolve
from multiprocessing import Pool, cpu_count
# from Funcs.ode_multip import solve_ode
from Funcs.ode_solv import ode_solv
from Funcs.plots_updating import plot_concentration_box_timeseries
import csv

def model_box(numLoop,  key_spe_for_plot, dt, rowvals, colptrs,plot_spec, formula, c, rrc,modelparams):
    comp_plot_indices = [modelparams.comp_namelist.index(plot_spec[i]) for i in range(len(plot_spec))]
    key_index = modelparams.comp_namelist.index(key_spe_for_plot)
    # Indices of species held at constant concentration
    _const_set = set(modelparams.const_comp)
    const_indices = np.array(
        [i for i, name in enumerate(modelparams.comp_namelist) if name in _const_set],
        dtype=int)

    delta_c_final = []
    tim_1_final = []
    c_final = []
    for k in range(numLoop):
        tim = round((k + 1) * modelparams.timesteps *dt,3)    
        old = c[key_index].copy()
        c= ode_solv(c, dt*modelparams.timesteps, rrc,modelparams.rindx_g, modelparams.pindx_g, modelparams.rstoi_g, modelparams.pstoi_g, modelparams.nreac_g, modelparams.nprod_g,
                    modelparams.jac_stoi_g, modelparams.njac_g, modelparams.jac_den_indx_g, modelparams.jac_indx_g,
                    modelparams.y_arr_g, modelparams.y_rind_g, modelparams.uni_y_rind_g, modelparams.y_pind_g, modelparams.uni_y_pind_g,
                    modelparams.reac_col_g, modelparams.prod_col_g, modelparams.rstoi_flat_g, modelparams.pstoi_flat_g, modelparams.rr_arr_g,
                    modelparams.rr_arr_p_g, rowvals, colptrs, modelparams.comp_num,modelparams.dil_fac_now,modelparams.wall_loss,
                    modelparams.jac_flat_rr_idx, modelparams.jac_flat_stoi, modelparams.jac_flat_den_indx, modelparams.jac_flat_data_indx,
                    const_comp_indices=const_indices)

        ######## set the constant concentrations for const_comp 
        for i in modelparams.const_comp:  
            try:
                c[modelparams.comp_namelist.index(i)] = getattr(modelparams, i+'conc')[0][0]
            except:
                pass
        
        ######## restore the constant concentration again        
        new = c[key_index]     
        delta_c = np.linalg.norm(new - old) / np.linalg.norm(old)
        tim_1 = round((k + 1) * modelparams.timesteps * dt, 3)
        delta_c_final.append(delta_c)
        tim_1_final.append(tim_1)
        c_final.append(c.copy())

        if k % modelparams.num_plot == 0:
            plot_concentration_box_timeseries(c_final, tim_1_final,plot_spec, formula,  comp_plot_indices)
            print('t = ' + str(tim) + ' ' + str(key_spe_for_plot) + " difference ratio: {:.4f} ".format(delta_c)+str(key_spe_for_plot)+" conc: {:.2E}".format(new))
        if k > 50 and delta_c < modelparams.actual_acc:
            break 

    delta_path = f"{modelparams.export_file_folder}{modelparams.file_name[:-4]}R{modelparams.Rgrid}L{modelparams.Zgrid}inter_acc{modelparams.actual_acc:.0e}dt{modelparams.dt:.0e}timstep{modelparams.timesteps:.0e}itx{modelparams.Itx:.0e}{modelparams.model_mode}.csv"
    with open(delta_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(modelparams.comp_namelist + ["time", "delta_c"])

        for comp_list, t, delta_val in zip(c_final, tim_1_final, delta_c_final):
            row = list(comp_list) + [t, delta_val]
            writer.writerow(row)

    return c_final
