# SA (sulfuric acid) calibration example — MCM mechanism
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

    Rl = 0.78, # ID for the tube with UV lights on
    Ll = 2.7,  # length for the tube with UV lights on 

    R1 = 0.78, # ID for the 1st tube
    L1 = 26,  # length for the 1st tube
    R2 = 1.2,
    L2 = 68,
    Itx = 1.8e-10, #1.25e-10,  # IT product at Qx; if unavailable, it needs to be calculated and conduct it product experiment see 10.5194/amt-4-437-2011
    Qx = 20,  # Qx where the IT product was calculated
    sampleflow = 22.4, # slpm

    outflowLocation = 'before',  # Outflow tube location: 'before' or 'after' injecting air, water, and SO2
    fullOrSimpleModel = 'full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (slower)

    ISOPratio = 5e-6,  # Isoprene ratio of the gas bottle (in ppm)
    SO2ratio = 5e-3,  # SO2 ratio of the gas bottle (in ppm)
    O2ratio = 0.209,  # O2 ratio in synthetic air
    ### if you have the OH source as continuous, you need to set the grid for this part
    Zgridl = 6, 
    Rgridl = 30,  
    ### this grid is total grid for the latter simulation
    Zgrid1 = 10,  # Number of grid points along tube length (commonly 80)
    Rgrid1 = 30,  # Number of grid points in radius direction (commonly 40)
    ### this grid is total grid for the latter simulation
    Zgrid2 = 10,  # Number of grid points along tube length (commonly 80)
    Rgrid2 = 30,  # Number of grid points in radius direction (commonly 40)

    dt = 2e-3, # Differential time interval (s) default is 2e-3
    timesteps = 8000,  # Number of time steps for integration for diffusion and convection part 
    inter_acc = 1e-4, # the difference of concentration between the previous interation and new one 1e-5 for flowtube1, higher for flowtube2 , but needs to div 1e3
    fix_timstep = 2, # the time steps from breaking the 1st simulation to the next twotubes. 
    model_mode= 'flowtube2',  # 'kinetic', 'box', or 'flowtube' (default: flowtube) 'kinetic': diffusion and conversion but without chemistry; flowtube1: chemisty equations are less than 50, flowtube2: chemistry equations are more than 50 and if you have continous OH source; # box: without diffusion and conversion 

    num_plot = 40, # Plot every 40 integrations
    # Species with user-defined diffusion values; otherwise, diffusion is calculated automatically
    Diff_setname = ['OH', 'HO2', 'SO3', 'H2SO4', 'H', 'O'],
    Diff_set = [0.215, 0.141, 0.126, 0.088, 0.20, 0.20],

    # Chemical mechanism files
    sch_name='mcm_export.fac',
    tsv_file='mcm_export_species.tsv', # File containing molecular formulas and SMILES
    folder_mechaism = 'SA',
    flag_mech = '1', # 1: MCM, 0: other mechanism

    # Constant concentration species
    const_comp=['SO2','O2','H2O'],

    # Species initialized at the first grid point
    Init_comp = ['OH','HO2'],
    Init_set  = 'off',
    # Key species for stopping criteria
    key_spe_for_plot = 'H2SO4',

    # Species to be plotted
    plot_spec=['OH','SO3','H2SO4','HO2','HSO3','H2O2','SO2','H2O','O2'],

    # Output filename (also used as input CSV from Input_files/)
    file_name='OH_it_production.csv',

    # Box model parameters
    dil_fac_now = 0,  # Dilution factor (used in box model)
    wall_loss_set = 0,  # Wall loss factor (used in box model)

    # Simulation time settings
    tot_time = 10,  # Total integration time (s) 10
    save_step = 0.5,  # Time interval for saving results (s) 0.5

    # Skip chemistry calculations (set to 1 to skip)
    pars_skip = 0,  
)

'''''''''
Calculate the input concentrations based on the parameters
'''''''''
calculate_concentrations(modelparams)
# Add CO at 100 ppb as constant background species
CO_conc = 100e-9 * modelparams.p / (1.3806488e-23 * modelparams.TEMP * 1e6)
modelparams.COconc = np.full_like(modelparams.O2conc, CO_conc)
modelparams.const_comp = list(modelparams.const_comp) + ['CO']
modelparams.const_comp_conc = np.transpose([getattr(modelparams, i + 'conc') for i in modelparams.const_comp])
'''''''''
prepare the input concentration
'''''''''
if modelparams.model_mode == 'kinetic':  # This mode runs the code in kinetic mode in which chemistry does not exist
    # define the H2SO4 concentration as 1e8 for convenience. !!!! This needs improvement
    OHconc = np.full_like(modelparams.OHconc, 1e8, dtype=np.float32)  # Set H2SO4 as 1e8
    modelparams.Init_comp = ['OH', 'HO2']  # here you need to change the initial compounds that you already set in paras
    modelparams.Init_comp_conc = np.column_stack([OHconc] * len(modelparams.Init_comp))
'''''''''
Run Flowtube Model
'''''''''
if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")  # Explicitly set the method for safety
    num_stage = modelparams.OHconc
    # Start timing execution
    start_time = time.perf_counter()

    Run_flowtube(modelparams,  num_stage)

    # End timing execution
    end_time = time.perf_counter()  # End timer
    total_time = (end_time - start_time)/60 # Compute elapsed time
    # Ensure print happens only in the main process
    if multiprocessing.current_process().name == "MainProcess":
        print(f"Total execution time: {total_time:.2f} mins")
    modelparams.total_time = total_time  # Store total time in the 
