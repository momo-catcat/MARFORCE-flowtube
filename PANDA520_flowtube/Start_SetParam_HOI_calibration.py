# =====================================================================================
# Example 2: HOI calibration with a POINT OH source
# =====================================================================================
# OH and HO2 from H2O photolysis at the inlet (as in Example 1) react with I2:
#     I2 + OH -> HOI + I   (+ iodine and HOx side reactions, input_mechanism/HOI/HOI_cali_chem.txt).
# The measured I2 of every stage (column I2conc of the input CSV) is held constant along the tube.
# The tube consists of two sections of different radius (R1/L1 then R2/L2) carrying the same flow;
# the model simulates them as two tubes and hands the outlet profile of the 1st to the inlet of the 2nd.
#
# Input: Input_files/HOI_calibration.csv (2 stages).
# Run from PANDA520_flowtube/:   python Start_SetParam_HOI_calibration.py
# =====================================================================================
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

        self.Rgrid = int(self.Rgrid1)  # Number of radial grid points used by the solver
        # Convert all lists, tuples, or NumPy arrays of numbers to np.float32 arrays
        for key, value in self.__dict__.items():
            if isinstance(value, (list, tuple, np.ndarray)):  # Check if value is a list, tuple, or array
                if all(isinstance(i, (int, float)) for i in value):  # Ensure all elements are numeric
                    setattr(self, key, np.array(value, dtype=np.float32))  # Convert to np.float32 array

modelparams = SimulationParams(
    # --- Experimental conditions ---------------------------------------------------
    p = 101000,                 # Pressure (Pa)
    TEMP = 298,                 # Temperature (K)

    # --- Tube geometry. No Rl/Ll -> point OH source ----------------------------------
    R1 = 0.78,                  # Inner RADIUS of the 1st section (cm)
    L1 = 41,                    # Length of the 1st section (cm)
    R2 = 1.04,                  # Inner RADIUS of the 2nd section (cm)
    L2 = 58.5,                  # Length of the 2nd section (cm)

    # --- OH source (point) -------------------------------------------------------------
    Itx = 4.84e10,              # Lamp It product (photons cm-2) measured at flow Qx; see doi:10.5194/amt-4-437-2011
    Qx = 20,                    # Flow (slpm) at which Itx was measured

    # --- Flows ---------------------------------------------------------------------------
    sampleflow = 22.5,          # CIMS inlet flow (slpm)
    outflowLocation = 'before', # Exhaust 'before' or 'after' the gas injection (see Example 1)
    O2ratio = 0.209,            # O2 fraction of the synthetic air in the O2flow column

    # --- Grid: tube 1 (Rgrid1 x Zgrid1 over L1) and tube 2 (Rgrid2 x Zgrid2 over L2) -----
    Rgrid1 = 80,                # Radial grid points across the diameter, tube 1 (80 = original resolution)
    Zgrid1 = 20,                # Axial grid points, tube 1
    Rgrid2 = 80,                # Radial grid points, tube 2 (keep equal to Rgrid1)
    Zgrid2 = 20,                # Axial grid points, tube 2

    # --- Numerics ----------------------------------------------------------------------
    model_mode = 'flowtube1',   # Chemistry + diffusion + advection solved together in every time step
    dt = 1e-4,                  # Time step (s). The HOI chemistry is slow (fastest loss: I2 + OH, ~1 s-1),
                                # so dt is set by the stability of diffusion on the fine R80 grid
                                # (dr = 0.02 cm): dt = 2e-4 already gives NaN here, 1e-4 is stable.
    timesteps = 5000,           # Time steps per iteration; one iteration = dt * timesteps = 0.5 s
    inter_acc = 1e-4,           # Converged when HOI at the outlet changes by < 0.01 % between iterations
    fix_timstep = 2,            # Minimum number of iterations before convergence is accepted
    num_plot = 2,               # Print progress / update the plots every 2 iterations

    # --- Diffusion coefficients (cm2 s-1); the others are estimated from their formula
    Diff_setname = ['OH', 'HO2'],
    Diff_set     = [0.215, 0.141],

    # --- Chemical mechanism (input_mechanism/HOI/) --------------------------------------
    folder_mechaism = 'HOI',                    # Sub-folder of input_mechanism/
    sch_name = 'HOI_cali_chem.txt',             # Mechanism file ('% rate : reaction ;' format)
    tsv_file = 'chemical_species_custom.xml',   # Non-MCM mechanism: species/SMILES are read from
                                                # chemical_species_custom.xml in the same folder
    flag_mech = '0',                            # '1' = MCM export, '0' = other mechanism

    # --- Species treatment -------------------------------------------------------------
    const_comp = ['I2', 'O2', 'H2O'],   # Held constant everywhere; I2 from the I2conc column
    Init_comp = ['OH', 'HO2'],          # Fixed at the tube inlet (point OH source) ...
    Init_set  = 'on',                   # ... 'on' = keep them fixed at the inlet
    key_spe_for_plot = 'HOI',           # Species used for the convergence test and progress output
    plot_spec = ['OH', 'HOI', 'HO2', 'I', 'I2'],  # Species shown in the plots

    final_output_method = 'mean',       # Outlet concentration: 'mean' = area average over the cross-section
                                        # (default), 'weighted' = flow-weighted (what leaves the tube per
                                        # unit time, as in the Matlab model); see README section 7
    # --- Input / output ----------------------------------------------------------------
    file_name = 'HOI_calibration.csv',  # Input CSV in Input_files/; results in Export_files/HOI_calibration/

    # --- Box-model and internal settings (not used for this flow-tube run; keep as is) -
    dil_fac_now = 0,
    wall_loss_set = 0,
    tot_time = 10,
    save_step = 0.5,
    pars_skip = 0,              # 0 = read and parse the mechanism file (normal)
)

'''''''''
Calculate the concentrations of every stage from the flows in the input CSV
'''''''''
calculate_concentrations(modelparams)

'''''''''
Run the model
'''''''''
if __name__ == "__main__":
    multiprocessing.set_start_method("spawn")   # required for the parallel solver on macOS/Windows
    num_stage = modelparams.OHconc              # runs every stage (CSV row) with OH > 0
    start_time = time.perf_counter()

    Run_flowtube(modelparams, num_stage)

    total_time = (time.perf_counter() - start_time) / 60
    if multiprocessing.current_process().name == "MainProcess":
        print(f"Total execution time: {total_time:.2f} mins")
