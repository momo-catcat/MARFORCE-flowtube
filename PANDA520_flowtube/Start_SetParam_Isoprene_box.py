# =====================================================================================
# Example 3: isoprene chemistry in a 0-D BOX MODEL (no transport)
# =====================================================================================
# Well-mixed chemistry only (Wennberg reduced isoprene mechanism, same file as Example 4).
# OH, CO, ISOP1OH2OOH, O2 and H2O are held constant; first-order dilution and species-specific
# wall losses mimic a chamber (values from a FOAM box-model set-up). The chemistry is integrated
# in chunks of dt * timesteps until the key species reaches steady state.
# Box mode uses only the FIRST row of the input CSV (one stage).
# SUN photolysis is off here (sun_a not set -> SUN = 0).
#
# Input: Input_files/Isoprene_box.csv.
# Run from PANDA520_flowtube/:   python Start_SetParam_Isoprene_box.py
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
    p = 101300,                 # Pressure (Pa)
    TEMP = 293.15,              # Temperature (K)

    # --- Lamp / geometry -----------------------------------------------------------------
    # Not used for transport in box mode, but Rl > 0 switches on the lamp reaction of the mechanism
    # (itProd: H2O = OH + H, rate coefficient Itx), which adds H -> HO2. OH itself is held constant.
    Rl = 0.78,                  # Inner RADIUS of the illuminated section (cm)
    Ll = 2.7,                   # Length of the illuminated section (cm)
    R1 = 0.78,                  # Inner RADIUS of the 1st tube (cm)
    L1 = 26,                    # Length of the 1st tube (cm)
    R2 = 1.2,                   # Inner RADIUS of the 2nd tube (cm)
    L2 = 68,                    # Length of the 2nd tube (cm)
    Itx = 1.8e-10,              # First-order rate coefficient of H2O -> OH + H (s-1)
    Qx = 20,                    # Flow (slpm) at which Itx was determined
    sampleflow = 22.4,          # CIMS inlet flow (slpm)
    outflowLocation = 'before', # Exhaust 'before' or 'after' the gas injection
    O2ratio = 0.209,            # O2 fraction of the synthetic air

    # --- Grid (required internally, not used by the box model) ---------------------------
    Zgridl = 6,
    Rgridl = 14,
    Zgrid1 = 10,
    Rgrid1 = 14,
    Zgrid2 = 10,
    Rgrid2 = 14,

    # --- Numerics ----------------------------------------------------------------------
    model_mode = 'box',         # 0-D box model: chemistry, dilution and wall loss only
    dt = 5e-3,                  # Together with timesteps: the chemistry is integrated (stiff solver) in
    timesteps = 8000,           # chunks of dt * timesteps = 40 s; dt itself has no stability limit here
    inter_acc = 1e-10,          # Box mode stops (after > 50 chunks) when the key species changes by less than
                                # inter_acc * timesteps / 1e4 (= 8e-11) between chunks, i.e. at steady state
    fix_timstep = 2,            # (flow-tube setting, not used in box mode)
    num_plot = 40,              # Print progress / update the time-series plot every 40 chunks

    # --- Diffusion coefficients (not used in box mode) -----------------------------------
    Diff_setname = ['OH', 'HO2'],
    Diff_set     = [0.215, 0.141],

    # --- Chemical mechanism (input_mechanism/Isoprene/Wennberg/) ------------------------
    folder_mechaism = 'Isoprene/Wennberg',
    sch_name = 'isoprene_reduced_plus_v5_SUNfinal.eqn',  # Wennberg reduced mechanism (KPP format)
    tsv_file = 'SMILE_formula.csv',     # Species SMILES table (species XML is chemical_species_custom.xml)
    flag_mech = '0',                    # '1' = MCM export, '0' = other mechanism

    # --- Species treatment -------------------------------------------------------------
    # Held constant; concentrations from the *conc columns of the input CSV (OHconc, COconc,
    # ISOP1OH2OOHconc) and from the O2/H2O flows
    const_comp = ['O2', 'H2O', 'CO', 'ISOP1OH2OOH', 'OH'],
    Init_comp = [],                     # No inlet species in box mode
    Init_set  = 'off',
    key_spe_for_plot = 'IDHDP',         # Species used for the steady-state test and progress output
    plot_spec = ['ISOP1OH2OOH', 'H2O', 'OH', 'HO2', 'IDHDP', 'ICPDH', 'IDHPE', 'HAC', 'GLYC',
                 'IHPOO1', 'IHPOO2', 'IHPOO3'],

    final_output_method = 'mean',       # Outlet concentration: 'mean' = area average over the cross-section
                                        # (default), 'weighted' = flow-weighted (what leaves the tube per
                                        # unit time, as in the Matlab model); see README section 7
    # --- Input / output ----------------------------------------------------------------
    file_name = 'Isoprene_box.csv',     # Input CSV in Input_files/; results in Export_files/Isoprene_box/

    # --- Box-model losses ----------------------------------------------------------------
    dil_fac_now = 2.1e-4,               # First-order dilution rate (s-1), applied to all non-constant species
    # Alternatively compute the dilution from the chamber: k_dil = (total_flow / 60) / chamber_volume
    # chamber_volume = 26500,           # Chamber volume (L)
    # total_flow = 200,                 # Total flow (slpm)
    wall_loss_set = 0,                  # 0 = no general wall loss, 1 = estimate for all species, x = x s-1 for all
    wall_loss_custom = {                # Species-specific first-order wall loss (s-1), overrides wall_loss_set
        'IHOO4': 1.6e-3, 'IHOO1': 1.6e-3,
        'ISOP1OH2OOH': 1.6e-3, 'ISOP3OOH4OH': 1.6e-3,
        'ISOP1OH4OOH': 1.6e-3, 'ISOP1OOH4OH': 1.6e-3,
        'SA': 1.6e-3, 'HO2': 1.6e-3, 'OH': 1.6e-3,
        'IPNOONO2': 1.6e-3,
    },

    # --- Internal settings (keep as is) ------------------------------------------------
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
