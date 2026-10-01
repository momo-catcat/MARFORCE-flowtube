import numpy as np

def grid_para(modelparams):

    if (modelparams.OHsource == 'Continuous'):
        modelparams.Zgrid = int(modelparams.Zgrid1 + modelparams.Zgridl)  # 
        modelparams.dx = np.zeros([int(modelparams.Rgrid), int(modelparams.Zgrid), modelparams.comp_num], dtype=np.float32)
        # modelparams.dr = np.zeros([int(modelparams.Rgrid), int(modelparams.Zgrid), modelparams.comp_num], dtype=np.float32)
        if modelparams.R2 == 0:
            modelparams.R2 = modelparams.R1

        modelparams.dx[:, 0:modelparams.Zgridl, :] = (modelparams.Ll) / (modelparams.Zgridl-1)
        modelparams.dx[:, modelparams.Zgridl:, :] = (modelparams.L1) / (modelparams.Zgrid1- 1)
        modelparams.dr = 2 * modelparams.Rl / (modelparams.Rgridl - 1)
    
    elif (modelparams.OHsource == 'point') and (modelparams.Q1[0] != modelparams.Q2[0]):

        modelparams.Zgrid = int(modelparams.Zgrid1) 
        modelparams.dx = np.zeros([int(modelparams.Rgrid), int(modelparams.Zgrid), modelparams.comp_num], dtype=np.float32)
        # modelparams.dr = np.zeros([int(modelparams.Rgrid), int(modelparams.Zgrid), modelparams.comp_num], dtype=np.float32)

        modelparams.dx[:, :, :] = (modelparams.L1) / (modelparams.Zgrid1-1)
        modelparams.dr = 2 * modelparams.R1 / (modelparams.Rgrid1 - 1)

    elif (modelparams.Q1[0] == modelparams.Q2[0]):
        modelparams.Zgrid = int(modelparams.Zgrid1+ modelparams.Zgrid2)  #
        modelparams.dx = np.zeros([int(modelparams.Rgrid), int(modelparams.Zgrid), modelparams.comp_num], dtype=np.float32)
        # modelparams.dr = np.zeros([int(modelparams.Rgrid), int(modelparams.Zgrid), modelparams.comp_num], dtype=np.float32)

        modelparams.dx[:, :, :] = (modelparams.L1+modelparams.L2) / (modelparams.Zgrid-1)
        modelparams.dr = 2 * modelparams.R1 / (modelparams.Rgrid1 - 1)
    return modelparams
