# Example 3: isoprene chemistry in a 0-D BOX MODEL (no diffusion or flow)
#
# Same mechanism as Example 4. OH, CO, ISOP1OH2OOH, O2 and H2O are held constant; dilution and
# species-specific wall losses mimic a chamber (values taken from a FOAM box-model set-up).
# Box mode runs the first row of the input CSV only (one stage). SUN photolysis is off (sun_a not set).
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
    p = 101300, # Pressure, Pa (1013 mbar, matching FOAM)
    TEMP = 293.15,  # Temperature, K

    Rl = 0.78, # ID for the tube with UV lights on
    Ll = 2.7,  # length for the tube with UV lights on

    R1 = 0.78, # ID for the 1st tube
    L1 = 26,  # length for the 1st tube
    R2 = 1.2,
    L2 = 68,
    Itx = 1.8e-10,  # IT product at Qx; calibrated value
    Qx = 20,  # Qx where the IT product was calculated
    sampleflow = 22.4, # slpm

    outflowLocation = 'before',  # Outflow tube location: 'before' or 'after' injecting air, water, and SO2
    fullOrSimpleModel = 'full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (slower)

    ISOPratio = 5e-6,  # Isoprene ratio of the gas bottle (in ppm)
    SO2ratio = 5e-3,  # SO2 ratio of the gas bottle (in ppm)
    O2ratio = 0.209,  # O2 ratio in synthetic air
    ### Grid parameters (still needed internally, but box model ignores spatial dimensions)
    Zgridl = 6,
    Rgridl = 14,
    Zgrid1 = 10,
    Rgrid1 = 14,
    Zgrid2 = 10,
    Rgrid2 = 14,

    dt = 5e-3, # Differential time interval (s)
    timesteps = 8000,  # Number of time steps for integration
    inter_acc = 1e-10, # convergence threshold
    fix_timstep = 2,
    model_mode = 'box',  # Box model: no diffusion, no flow, well-mixed chemistry only

    num_plot = 40, # Plot every 40 integrations
    # Species with user-defined diffusion values; otherwise, diffusion is calculated automatically
    Diff_setname = ['OH', 'HO2', 'SO3', 'H2SO4'],
    Diff_set = [0.215, 0.141, 0.126, 0.088],

    # Chemical mechanism files
    sch_name = 'isoprene_reduced_plus_v5_SUNfinal.eqn',  # same mechanism as Example 4
    tsv_file = 'SMILE_formula.csv',
    folder_mechaism = 'Isoprene/Wennberg',
    flag_mech = '0', # 1: MCM, 0: other mechanism

    # Constant concentration species
    const_comp = ['O2', 'H2O', 'CO', 'ISOP1OH2OOH', 'OH'],

    # No Init_comp needed — all reactive species are in const_comp
    Init_comp = [],
    Init_set  = 'off',

    # Key species for stopping criteria
    key_spe_for_plot = 'IDHDP',

    # Species to be plotted
    plot_spec = ['ISOP1OH2OOH','H2O','OH','HO2','IDHDP','ICPDH','IDHPE','HAC','GLYC','IHPOO1','IHPOO2','IHPOO3'],

    # Output filename
    file_name = 'Isoprene_box.csv',

    # ── Box model parameters ──
    dil_fac_now = 2.1e-4,  # Dilution factor (s⁻¹) from FOAM
    wall_loss_set = 0,  # 0: no wall loss, 1: auto-calculate for all species
    wall_loss_custom = {
        'IHOO4': 1.6e-3, 'IHOO1': 1.6e-3,
        'ISOP1OH2OOH': 1.6e-3, 'ISOP3OOH4OH': 1.6e-3,
        'ISOP1OH4OOH': 1.6e-3, 'ISOP1OOH4OH': 1.6e-3,
        # IEPOXt/IEPOXc/IEPOXD wall loss removed to match FOAM effective behavior
        # (FOAM mechanism has wall loss reactions but output shows they are not applied)
        'SA': 1.6e-3, 'HO2': 1.6e-3, 'OH': 1.6e-3,
        'IPNOONO2': 1.6e-3,
    },  # Species-specific wall loss (s⁻¹) from FOAM

    # Auto-calculate dilution from chamber volume and flow (uncomment to use):
    # chamber_volume = 26500,  # Chamber volume in liters
    # total_flow = 200,        # Total flow in slpm
    # → k_dil = (total_flow / 60) / chamber_volume  [s⁻¹]

    # Simulation time settings
    tot_time = 10,  # Total integration time (s)
    save_step = 0.5,  # Time interval for saving results (s)

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
if modelparams.model_mode == 'kinetic':
    OHconc = np.full_like(modelparams.OHconc, 1e8, dtype=np.float32)
    modelparams.Init_comp = ['OH', 'HO2']
    modelparams.Init_comp_conc = np.column_stack([OHconc] * len(modelparams.Init_comp))

'''''''''
Run Box Model
'''''''''
if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")
    num_stage = modelparams.OHconc
    start_time = time.perf_counter()

    Run_flowtube(modelparams, num_stage)

    end_time = time.perf_counter()
    total_time = (end_time - start_time) / 60
    if multiprocessing.current_process().name == "MainProcess":
        print(f"Total execution time: {total_time:.2f} mins")
    modelparams.total_time = total_time
