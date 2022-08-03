import numpy as np
from matplotlib import pyplot as plt

from meanconc_cal import meanconc_cal_sim as meanconc_cal_sim
from odesolve3 import odesolve as odesolve
from odesolve3_Y import odesolve as odesolve_Y


def model_4(R2, R1, Rgrid, Zgrid, comp_num, L2, L1, numLoop, comp_namelist, key_spe_for_plot, dt, timesteps, Diff_vals,
            Rtot, const_comp_free, const_comp_conc_free, Q1, Q2, dydt_vst, rindx, nreac, rstoi, rate_values, const_comp, \
            u, plot_spec, formula, c, dr, dx):
    # %% first tube run
    for j in range(numLoop):
        c1 = c.copy()
        old = c1[:, -1, comp_namelist.index(key_spe_for_plot)]

        c = odesolve(timesteps, Zgrid, Rgrid, dt, Diff_vals, Rtot, dr, dx, Q1, c, comp_namelist, dydt_vst,
                     rindx,
                     nreac, rstoi, rate_values, const_comp, u)

        new = c[:, -1, comp_namelist.index(key_spe_for_plot)]

        if (j > 5) & (np.sum(new - old) / np.sum(old) < 1e-5):
            break
    # %% transfer the flow distribution for next run
    meanConc = meanconc_cal_sim(R1, Rgrid, comp_namelist, c)

    c = np.zeros([Rgrid, Zgrid, comp_num])
    for i in range(len(comp_namelist)):
        c[:int(Rgrid / 2), 0, i] = meanConc[i]
        c[int(Rgrid / 2):, 0, i] = 0
    for i in const_comp:
        c[:int(Rgrid / 2), :, comp_namelist.index(i)] = meanConc[comp_namelist.index(i)]

    # c[int(Rgrid / 2):, 0, :] = np.zeros([int(Rgrid / 2), comp_num])
    for i in const_comp_free:
        c[int(Rgrid / 2):, :, comp_namelist.index(i)] = const_comp_conc_free[const_comp_free.index(i)]

    # %% second tube run

    dx = L2 / (Zgrid - 1)
    dr = 2 * R2 / (Rgrid - 1)
    for j in range(numLoop):
        c1 = c.copy()
        old = c1[:, -1, comp_namelist.index(key_spe_for_plot)]

        c = odesolve_Y(timesteps, Zgrid, Rgrid, dt, Diff_vals, R2, dr, dx, Q2, c, comp_namelist, dydt_vst, rindx,
                        nreac, rstoi, rate_values, const_comp, u)

        tim = (j + 1) * timesteps * dt
        comp_plot_index = [comp_namelist.index(plot_spec[i]) for i in range(len(plot_spec))]

        new = c[:, -1, comp_namelist.index(key_spe_for_plot)]

        fig, axs = plt.subplots(2, 3, figsize=(8, 5), facecolor='w', edgecolor='k')
        fig.subplots_adjust(hspace=.5, wspace=.45)
        plt.style.use('default')
        plt.rcParams.update(
               {'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif', 'font.serif': 'Times New Roman'})

        axs = axs.ravel()

        for i in range(len(plot_spec)):
            axs[i].pcolor(np.linspace(0, L2 + L1, Zgrid), np.linspace(-R2, R2, Rgrid), c[:, :, comp_plot_index[i]],
                          shading='nearest', cmap='jet')

            axs[i].set_xlabel('L [cm]')
            axs[i].set_ylabel('R [cm]')
            axs[i].set_title(formula[i])

        fig.delaxes(axs[5])
        plt.gcf().text(0.7, 0.3, 'Time = ' + str(tim), fontsize=15)
        plt.draw()
        plt.pause(1)
        plt.close('all')
        print(['t = ' + str(tim) + str(key_spe_for_plot) + " difference: " + str(np.sum(new - old))])

        if (j > 5) & (np.sum(new - old) / np.sum(old) < 1e-5):
            break
    return c
