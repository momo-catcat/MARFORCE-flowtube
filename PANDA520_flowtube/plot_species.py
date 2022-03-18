def plot(c, L1, L2, R1, R2, species, Zgrid, Rgrid, comp_namelist):
    import numpy as np
    import matplotlib.pyplot as plt
    # sort chemistry   

    fig, axs = plt.subplots(1, 1, figsize=(5, 4), facecolor='w', edgecolor='k')
    plt.style.use('default')
    plt.rcParams.update(
        {'font.size': 13, 'font.weight': 'bold', 'font.family': 'serif', 'font.serif': 'Times New Roman'})
    axs.pcolor(np.linspace(0, L2 + L1, Zgrid), np.linspace(-R1, R1, Rgrid), c[:, :, comp_namelist.index(species)],
               shading='nearest')
    axs.pcolor(np.linspace(0, L2 + L1, Zgrid), np.linspace(-R1, R2, Rgrid), c[:, :, comp_namelist.index(species)],
               shading='nearest')
    # axs = axs.ravel()
    axs.set_xlabel('L [cm]')
    axs.set_ylabel('R [cm]')
    axs.set_title(species)

    plt.draw()
    plt.pause(1)
    print(species + ' last columns', c[:, -1, comp_namelist.index(species)])

    return (species)
