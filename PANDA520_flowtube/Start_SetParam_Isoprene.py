# Isoprene-OH oxidation example — Wennberg reduced mechanism
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
    Itx = 1.8e-10,  # IT product at Qx; sensitivity test: 10x lower than calibrated 1.8e-10
    Qx = 20,  # Qx where the IT product was calculated
    sampleflow = 22.4, # slpm

    outflowLocation = 'before',  # Outflow tube location: 'before' or 'after' injecting air, water, and SO2
    fullOrSimpleModel = 'full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (slower)

    # Add initial concentrations of species if you have 
    ISOPratio = 5e-6,  # Isoprene ratio of the gas bottle (in ppm)
    SO2ratio = 5e-3,  # SO2 ratio of the gas bottle (in ppm)
    O2ratio = 0.209,  # O2 ratio in synthetic air
    ### if you have the OH source as continuous, you need to set the grid for this part
    Zgridl = 6,
    Rgridl = 14,
    ### this grid is total grid for the latter simulation
    Zgrid1 = 10,  # Number of grid points along tube length (commonly 80)
    Rgrid1 = 14,  # Number of grid points in radius direction (commonly 40)
    ### this grid is total grid for the latter simulation
    Zgrid2 = 10,  # Number of grid points along tube length (commonly 80)
    Rgrid2 = 14,  # Number of grid points in radius direction (commonly 40)

    dt = 5e-4, # Differential time interval (s)
    timesteps = 8000,  # Number of time steps for integration
    inter_acc = 3e-4, # convergence threshold (matches converged R14 1.8e-10 reference run)
    fix_timstep = 2, # the time steps from breaking the 1st simulation to the next twotubes. 
    model_mode= 'flowtube2',  # 'kinetic', 'box', or 'flowtube' (default: flowtube) 'kinetic': diffusion and conversion but without chemistry; flowtube1: chemisty equations are less than 50, flowtube2: chemistry equations are more than 50 and if you have continous OH source; # box: without diffusion and conversion 

    num_plot = 40, # Plot every 40 integrations
    # Species with user-defined diffusion values; otherwise, diffusion is calculated automatically
    Diff_setname = ['OH', 'HO2', 'SO3', 'H2SO4'],
    Diff_set = [0.215, 0.141, 0.126, 0.088],  

    # Chemical mechanism files
    sch_name = 'isoprene_reduced_plus_v5.eqn',  # Chemical scheme file (stored in input_mechanism folder)
    tsv_file = 'SMILE_formula.csv',
    folder_mechaism = 'Isoprene/Wennberg',
    flag_mech = '0', # 0: non-MCM mechanism

    # Constant concentration species
    const_comp = ['O2', 'H2O'],

    # Species initialized at the first grid point
    Init_comp = ['ISOP1OH2OOH'],#['OH','HO2'],
    Init_set  = 'on',
    # Key species for stopping criteria
    key_spe_for_plot = 'IDHDP',  #'H2SO4',

    # Species to be plotted
    plot_spec = ['ISOP1OH2OOH','H2O','OH','HO2','IDHDP','ICPDH','IDHPE','HAC','GLYC','IHPOO1','IHPOO2','IHPOO3'],

    # Output filename (also used as input CSV from Input_files/)
    file_name = 'ISOP_OH_diff_OH.csv',

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
