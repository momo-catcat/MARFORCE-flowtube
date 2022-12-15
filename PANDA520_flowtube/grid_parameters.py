import numpy as np


def grid_para(Zgrid, Rgrid, R2, R1, L2, L1, comp_num):
    sp_line = int(Zgrid * L1 / (L2 + L1))
    if R2 == 0:
        R2 = R1
    dr = np.zeros([int(Rgrid), int(Zgrid), comp_num])
    dx = (L2 + L1) / (Zgrid - 1)
    Rtot = np.zeros([int(Rgrid), int(Zgrid), comp_num])
    Rtot[:, 0:sp_line, :] = R1
    Rtot[:, sp_line:, :] = R2
    dr[:, 0:sp_line, :] = 2 * R1 / (Rgrid - 1)
    dr[:, sp_line:, :] = 2 * R2 / (Rgrid - 1)
    return Rtot, dr, dx, sp_line
