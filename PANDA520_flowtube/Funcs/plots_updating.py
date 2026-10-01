import numpy as np
from matplotlib import pyplot as plt
import math

def plot_concentration_profiles(c, plot_spec, formula, L1, L2, Zgrid, Rgrid, tim, R2, comp_plot_indices):
    #### plots the cocurrent concentration profile

    #     # Reuse figure and axes to improve performance
    n_plots = len(plot_spec)
    n_rows = math.ceil(n_plots / 3)
    n_cols = min(3, n_plots)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(9, 5), facecolor='w', edgecolor='k')
    plt.cla()
    fig.subplots_adjust(hspace=.5, wspace=.35)
    plt.style.use('default')
    plt.rcParams.update(
        {'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif', 'font.serif': ['DejaVu Serif']})
    
    x = np.linspace(0, L1 + L2, Zgrid)
    y = np.linspace(-R2, R2, Rgrid)

    axs = axs.ravel()
    for i, (ax, spec_idx) in enumerate(zip(axs, comp_plot_indices)):
        cl = ax.pcolor(x,y, c[:, :, spec_idx],
                    shading='nearest', cmap='jet')
        ax.set_title(formula[i],fontsize = 12)
        clb = plt.colorbar(cl, ax = axs[i])
        clb.formatter.set_powerlimits((0, 0))
        clb.formatter.set_useMathText(True)

    fig.supxlabel('L [cm]', fontsize=14)
    fig.supylabel('R [cm]', fontsize=14)
    fig.suptitle(f'Time = {tim:.3f}')
    plt.draw()
    plt.pause(0.01) # Non-blocking plot
    plt.close(fig)

def  plot_concentration_profiles_3(c_1st, c_2nd, plot_spec, formula, tim, comp_plot_indices, modelparams):

    # 1. Generate x (axial coordinate)
    Zl = int(modelparams.Zgridl)
    Z1 = int(modelparams.Zgrid1)
    Z2 = int(modelparams.Zgrid2)

    x_s = np.linspace(0, modelparams.Ll, Zl, endpoint=False)
    x1 = np.linspace(modelparams.Ll, modelparams.Ll + modelparams.L1, Z1, endpoint=False)
    x2 = np.linspace(modelparams.Ll + modelparams.L1, modelparams.Ll + modelparams.L1 + modelparams.L2, Z2, endpoint=False)
    x = np.concatenate((x_s, x1, x2))  # Total x: Zl + Z1 + Z2
    nr = c_1st.shape[0]  

    y = np.linspace(0, modelparams.Rl, nr).reshape(-1, 1)  # Uniform radial grid for all regions

    # Repeat radial coordinates along x
    y = np.repeat(y, Zl + Z1 + Z2, axis=1)  # shape: (nr, total_z)

    c_combined = np.concatenate((c_1st, c_2nd,), axis=1)
    # 3. Setup plots
    n_plots = len(plot_spec)
    n_rows = math.ceil(n_plots / 3)
    n_cols = min(3, n_plots)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(9, 5))
    fig.subplots_adjust(hspace=.5, wspace=.35)
    plt.style.use('default')
    plt.rcParams.update({'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif'})

    axs = axs.ravel()
    for i, (ax, spec_idx) in enumerate(zip(axs, comp_plot_indices)):
        cl = ax.pcolormesh(x, y, c_combined[:, :, spec_idx], shading='auto', cmap='jet')
        ax.set_title(formula[i], fontsize=12)
        clb = plt.colorbar(cl, ax=ax)
        clb.formatter.set_powerlimits((0, 0))
        clb.formatter.set_useMathText(True)
        clb.update_ticks()

    fig.supxlabel('Axial Length [cm]', fontsize=14)
    fig.supylabel('Radial Distance [cm]', fontsize=14)
    fig.suptitle(f'Time = {tim:.3f}')
    plt.draw()
    plt.pause(0.01)
    plt.close(fig)

def  plot_concentration_profiles_different_dx(c, plot_spec, formula, tim, comp_plot_indices, modelparams):

    # 1. Generate x (axial coordinate)
    Zl = int(modelparams.Zgridl)
    Z1 = int(modelparams.Zgrid1)

    x_s = np.linspace(0, modelparams.Ll, Zl, endpoint=False)
    x1 = np.linspace(modelparams.Ll, modelparams.Ll + modelparams.L1, Z1, endpoint=False)
    x = np.concatenate((x_s, x1))  # Total x: Zl + Z1 + Z2
    nr = c.shape[0]  

    y = np.linspace(0, modelparams.Rl, nr).reshape(-1, 1)  # Uniform radial grid for all regions

    # Repeat radial coordinates along x
    y = np.repeat(y, Zl + Z1, axis=1)  # shape: (nr, total_z)

    # 3. Setup plots
    n_plots = len(plot_spec)
    n_rows = math.ceil(n_plots / 3)
    n_cols = min(3, n_plots)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(9, 5))
    fig.subplots_adjust(hspace=.5, wspace=.35)
    plt.style.use('default')
    plt.rcParams.update({'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif'})

    axs = axs.ravel()
    for i, (ax, spec_idx) in enumerate(zip(axs, comp_plot_indices)):
        cl = ax.pcolormesh(x, y, c[:, :, spec_idx], shading='auto', cmap='jet')
        ax.set_title(formula[i], fontsize=12)
        clb = plt.colorbar(cl, ax=ax)
        clb.formatter.set_powerlimits((0, 0))
        clb.formatter.set_useMathText(True)
        clb.update_ticks()

    fig.supxlabel('Axial Length [cm]', fontsize=14)
    fig.supylabel('Radial Distance [cm]', fontsize=14)
    fig.suptitle(f'Time = {tim:.3f}')
    plt.draw()
    plt.pause(0.01)
    plt.close(fig)

def plot_concentration_box_timeseries(c_final, tim_1_final,plot_spec, formula,  comp_plot_indices):
    #### plots the cocurrent concentration profile

    #     # Reuse figure and axes to improve performance
    n_plots = len(plot_spec)
    n_rows = math.ceil(n_plots / 3)
    n_cols = min(3, n_plots)
    c_arr = np.array(c_final, dtype=float)
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(9, 5), facecolor='w', edgecolor='k')
    plt.cla()
    fig.subplots_adjust(hspace=.5, wspace=.35)
    plt.style.use('default')
    plt.rcParams.update(
        {'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif', 'font.serif': ['DejaVu Serif']})


    axs = axs.ravel()
    for i, (ax, spec_idx) in enumerate(zip(axs, comp_plot_indices)):
        y = c_arr[:, spec_idx]
        cl = ax.plot(tim_1_final, y, linewidth = 2)
        ax.set_title(formula[i],fontsize = 12)
        # clb = plt.colorbar(cl, ax = axs[i])
        # clb.formatter.set_powerlimits((0, 0))
        # clb.formatter.set_useMathText(True)

    fig.supxlabel('Time s', fontsize=14)
    fig.supylabel('Conc', fontsize=14)
    plt.draw()
    plt.pause(0.01) # Non-blocking plot
    plt.close(fig)

def plot_concentration_profiles_test(c_1st, plot_spec, formula, tim, comp_plot_indices, modelparams,ne):

    # 1. Generate x (axial coordinate)
    Zl = int(modelparams.Zgridl)
    Z1 = int(modelparams.Zgrid1)

    x_s = np.linspace(0, modelparams.Ll, Zl, endpoint=False)
    x1 = np.linspace(modelparams.Ll, modelparams.Ll + modelparams.L1, Z1, endpoint=False)
    x = np.concatenate((x_s, x1))  # Total x: Zl + Z1 + Z2
    nr = c_1st.shape[0]  

    y = np.linspace(0, modelparams.Rl, nr).reshape(-1, 1)  # Uniform radial grid for all regions

    # Repeat radial coordinates along x
    y = np.repeat(y, Zl + Z1 , axis=1)  # shape: (nr, total_z)

    # 3. Setup plots
    n_plots = len(plot_spec)
    n_rows = math.ceil(n_plots / 3)
    n_cols = min(3, n_plots)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(9, 5))
    fig.subplots_adjust(hspace=.5, wspace=.35)
    plt.style.use('default')
    plt.rcParams.update({'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif'})

    axs = axs.ravel()
    for i, (ax, spec_idx) in enumerate(zip(axs, comp_plot_indices)):
        cl = ax.pcolormesh(x, y, c_1st[:, :, spec_idx], shading='auto', cmap='jet')
        ax.set_title(formula[i], fontsize=12)
        clb = plt.colorbar(cl, ax=ax)
        clb.formatter.set_powerlimits((0, 0))
        clb.formatter.set_useMathText(True)
        clb.update_ticks()
    fig.suptitle(ne, fontsize=16)
    fig.supxlabel('Axial Length [cm]', fontsize=14)
    fig.supylabel('Radial Distance [cm]', fontsize=14)
    # fig.suptitle(f'Time = {tim:.3f}')
    plt.draw()
    plt.pause(0.01)
    plt.close(fig)

def plot_concentration_profiles_test1(c, plot_spec, formula, L1, L2, Zgrid, Rgrid, tim, R2, comp_plot_indices,be):
    #### plots the cocurrent concentration profile

    #     # Reuse figure and axes to improve performance
    n_plots = len(plot_spec)
    n_rows = math.ceil(n_plots / 3)
    n_cols = min(3, n_plots)

    fig, axs = plt.subplots(n_rows, n_cols, figsize=(9, 5), facecolor='w', edgecolor='k')
    plt.cla()
    fig.subplots_adjust(hspace=.5, wspace=.35)
    plt.style.use('default')
    plt.rcParams.update(
        {'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif', 'font.serif': ['DejaVu Serif']})
    
    x = np.linspace(0, L1 + L2, Zgrid)
    y = np.linspace(-R2, R2, Rgrid)

    axs = axs.ravel()
    for i, (ax, spec_idx) in enumerate(zip(axs, comp_plot_indices)):
        cl = ax.pcolor(x,y, c[:, :, spec_idx],
                    shading='nearest', cmap='jet')
        ax.set_title(formula[i],fontsize = 12)
        clb = plt.colorbar(cl, ax = axs[i])
        clb.formatter.set_powerlimits((0, 0))
        clb.formatter.set_useMathText(True)

    fig.supxlabel('L [cm]', fontsize=14)
    fig.supylabel('R [cm]', fontsize=14)
    fig.suptitle(f'Time = {tim:.3f}')
    fig.suptitle(be)
    plt.draw()
    plt.pause(0.01) # Non-blocking plot
    plt.close(fig)



