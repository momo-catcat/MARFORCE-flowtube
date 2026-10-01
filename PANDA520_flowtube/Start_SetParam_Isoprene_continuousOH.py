# Example 4: isoprene hydroxy-hydroperoxide (ISOP1OH2OOH) + OH — CONTINUOUS OH source, two tubes
#
# OH is produced all along the illuminated section (radius Rl, length Ll) by the mechanism reaction
#   % itProd(modelparams) : H2O = OH + H ;   with itProd = Itx (first-order rate coefficient, s-1)
# To use a CONTINUOUS source: set Rl / Ll / Rgridl / Zgridl and Init_set = "off" for OH (OH is not fixed at the inlet).
# ISOP1OH2OOH (column ISOP1OH2OOHconc of the input CSV) is fixed at the inlet (Init_set = "on").
#
# SUN photolysis: the mechanism contains reactions written as SUN*k (e.g. H2O2 = OH + OH : SUN*5e-5);
# SUN takes the value of sun_a below. Set sun_a = 0 to switch these photolysis reactions off.
# A constant 100 ppb CO background is added after calculate_concentrations().
import numpy as np
from Funcs.Calcu_by_flow import calculate_concentrations
import time
import multiprocessing
from Funcs.Run_flowtube import Run_flowtube
'''''''''
set parameters
'''''''''
class SimulationParams:
    """Class for storing all model parameters."""
    def __init__(self, **kwargs):
        """Dynamically load all parameters into self."""
        for key, value in kwargs.items():
            setattr(self, key, value)

        # Convert all numerical scalar parameters to np.float32
        for key, value in self.__dict__.items():
            if isinstance(value, (int, float)):  # Check if the value is a scalar number
                if not key.startswith("Zgrid") and not key.startswith("Rgrid") and key != "timesteps":
                    setattr(self, key, np.float32(value))

        for grid_key in ["Zgrid1", "Rgrid1", "Zgrid2", "Rgrid2", "Zgridl", "Rgridl",
            "Zgrid_p", "Zgrid_s", "Rgrid_p", "timesteps"]:
            if hasattr(self, grid_key):
                setattr(self, grid_key, int(getattr(self, grid_key)))

        self.Rgrid = int(self.Rgrid1)  # Total grid points along tube length
        # Convert all lists, tuples, or NumPy arrays of numbers to np.float32 arrays
        for key, value in self.__dict__.items():
            if isinstance(value, (list, tuple, np.ndarray)):  # Check if value is a list, tuple, or array
                if all(isinstance(i, (int, float)) for i in value):  # Ensure all elements are numeric
                    setattr(self, key, np.array(value, dtype=np.float32))  # Convert to np.float32 array

# Initialize model parameters
modelparams = SimulationParams(
    p = 101000, # Pressure, Pa
    TEMP = 293.15,  # Temperature, K

    # Illuminated section -> continuous OH source
    Rl = 0.78, # Inner radius of the illuminated section (cm)
    Ll = 2.7,  # Length of the illuminated section (cm)

    R1 = 0.78, # Inner radius of the 1st tube (cm)
    L1 = 26,  # Length of the 1st tube (cm)
    R2 = 1.2,  # Inner radius of the 2nd tube (cm)
    L2 = 68,   # Length of the 2nd tube (cm)
    Itx = 1.9e-10,  # Continuous source: first-order rate coefficient of H2O -> OH + H in the lamp section (s-1)
    sun_a = 2e5,  # Value of SUN in the mechanism (photolysis scaling); 0 = photolysis off
    Qx = 20,  # Flow (slpm) at which Itx was determined
    sampleflow = 22.4, # Inlet flow of the CIMS (slpm)

    outflowLocation = 'before',  # Outflow tube location: 'before' or 'after' injecting air, water, and isoprene
    fullOrSimpleModel = 'full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (slower)

    O2ratio = 0.209,  # O2 fraction in synthetic air

    # Grid: illuminated section (l), 1st tube (1), 2nd tube (2)
    Zgridl = 6,
    Rgridl = 14,
    Zgrid1 = 10,
    Rgrid1 = 14,
    Zgrid2 = 10,
    Rgrid2 = 14,

    dt = 2e-3, # Transport time step (s); reduce if you see NaN (finer grids need smaller dt)
    timesteps = 250,  # Transport steps per iteration; chemistry is integrated over dt * timesteps = 0.5 s
    inter_acc = 1e-3, # Convergence: relative change of key species between iterations
    fix_timstep = 2, # Minimum number of iterations
    model_mode = 'flowtube2',  # Large mechanism / continuous source: parallel chemistry solver

    num_plot = 5, # Plot/print every 5 iterations
    # Species with user-defined diffusion coefficients (cm2 s-1)
    Diff_setname = ['OH', 'HO2'],
    Diff_set = [0.215, 0.141],

    # Chemical mechanism (Wennberg reduced isoprene mechanism, KPP format)
    sch_name = 'isoprene_reduced_plus_v5_SUNfinal.eqn',
    tsv_file = 'SMILE_formula.csv',
    folder_mechaism = 'Isoprene/Wennberg',
    flag_mech = '0', # 1: MCM, 0: other mechanism

    # Constant concentration species (CO is added below)
    const_comp = ['O2', 'H2O'],

    # ISOP1OH2OOH fixed at the inlet
    Init_comp = ['ISOP1OH2OOH'],
    Init_set  = 'on',

    # Key species for the convergence criterion
    key_spe_for_plot = 'IDHDP',

    # Species to be plotted
    plot_spec = ['ISOP1OH2OOH', 'OH', 'HO2', 'IDHDP', 'ICPDH', 'IDHPE', 'HAC', 'GLYC', 'IHPOO1', 'IHPOO2', 'IHPOO3'],

    # Input CSV in Input_files/ (one row per stage); also the output folder name
    file_name = 'Isoprene_continuousOH.csv',

    # Box model parameters (not used in flow-tube mode)
    dil_fac_now = 0,
    wall_loss_set = 0,
    tot_time = 10,
    save_step = 0.5,

    pars_skip = 0,  # 0: parse the mechanism
)

'''''''''
Calculate the input concentrations from the flows
'''''''''
calculate_concentrations(modelparams)

# Add a constant 100 ppb CO background (any species can be added to const_comp this way)
CO_conc = 100e-9 * modelparams.p / (1.3806488e-23 * modelparams.TEMP * 1e6)
modelparams.COconc = np.full_like(modelparams.O2conc, CO_conc)
modelparams.const_comp = list(modelparams.const_comp) + ['CO']
modelparams.const_comp_conc = np.transpose([getattr(modelparams, i + 'conc') for i in modelparams.const_comp])

'''''''''
Run the model
'''''''''
if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    num_stage = modelparams.OHconc  # runs every stage (row) with OH > 0
    start_time = time.perf_counter()

    Run_flowtube(modelparams, num_stage)

    total_time = (time.perf_counter() - start_time) / 60
    if multiprocessing.current_process().name == "MainProcess":
        print(f"Total execution time: {total_time:.2f} mins")
