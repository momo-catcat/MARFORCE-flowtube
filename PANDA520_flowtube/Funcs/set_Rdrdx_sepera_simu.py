import numpy as np

def set_Rtot_dr_dx_for_spe_simulation(R2,L2,modelparams):
    dr2 = 2 * R2 / (modelparams.Rgrid2 - 1)
    dx2 = np.zeros([int(modelparams.Rgrid2), int(modelparams.Zgrid2), modelparams.comp_num], dtype=np.float32)
    dx2[:,:,:] = L2 / (modelparams.Zgrid2 - 1)

    return dr2, dx2


    # dx = L1 / (Zgrid - 1)
    # dr = np.zeros([int(Rgrid), int(Zgrid), modelparams.comp_num])
    # dr[:, :, :] = 2 * R1 / (Rgrid - 1)
    # Rtot = np.zeros([int(Rgrid), int(Zgrid), modelparams.comp_num])
    # Rtot[:, :, :] = R1