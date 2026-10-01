# =====================================================================================
# Example 4: ISOP1OH2OOH + OH in two tubes with a CONTINUOUS OH source and SUN photolysis
# =====================================================================================
# OH is produced all along the illuminated section (radius Rl, length Ll) by the mechanism reaction
#     % itProd(modelparams) : H2O = OH + H ;     with itProd = Itx (first-order rate coefficient, s-1),
# while ISOP1OH2OOH (column ISOP1OH2OOHconc of the input CSV) is fixed at the inlet.
# Tube 1 (R1, L1, flow Q1) is followed by tube 2 (R2, L2, flow Q2, which includes the flows added at the Y-piece).
#
# SUN photolysis: reactions written as SUN*k in the mechanism (e.g. SUN*5e-5 : H2O2 = OH + OH)
# use SUN = sun_a; set sun_a = 0 to switch them off.
# A constant 100 ppb CO background is added after calculate_concentrations().
#
# Continuous-source set-up = Rl/Ll/Rgridl/Zgridl set, Itx is a rate coefficient, OH not in
# Init_comp, model_mode = 'flowtube2'.  Restart files (checkpoints, warm start) are written for this
# case, so an interrupted run continues when the script is started again.
# Input: Input_files/Isoprene_continuousOH.csv (2 stages).
# Run from PANDA520_flowtube/:   python Start_SetParam_Isoprene_continuousOH.py
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
    TEMP = 293.15,              # Temperature (K)

    # --- Illuminated section -> continuous OH source -----------------------------------
    Rl = 0.78,                  # Inner RADIUS of the illuminated section (cm); Rl > 0 = continuous source
    Ll = 2.7,                   # Length of the illuminated section (cm)

    # --- Tubes ---------------------------------------------------------------------------
    R1 = 0.78,                  # Inner RADIUS of the 1st tube (cm)
    L1 = 26,                    # Length of the 1st tube (cm)
    R2 = 1.2,                   # Inner RADIUS of the 2nd tube (cm)
    L2 = 68,                    # Length of the 2nd tube (cm)

    # --- OH source (continuous) and photolysis -----------------------------------------
    Itx = 1.9e-10,              # First-order rate coefficient of H2O -> OH + H in the lamp section (s-1)
    Qx = 20,                    # Flow (slpm) at which Itx was determined
    sun_a = 2e5,                # Value of SUN in the mechanism (photolysis scaling); 0 = photolysis off

    # --- Flows ---------------------------------------------------------------------------
    sampleflow = 22.4,          # CIMS inlet flow (slpm)
    outflowLocation = 'before', # Exhaust 'before' or 'after' the gas injection
    O2ratio = 0.209,            # O2 fraction of the synthetic air in the O2flow columns

    # --- Grid: illuminated section (l), tube 1 (1), tube 2 (2) -------------------------
    Rgridl = 14,                # Radial points, illuminated section
    Zgridl = 6,                 # Axial points, illuminated section
    Rgrid1 = 14,                # Radial points, tube 1 (14 = fast example; publication runs used 60)
    Zgrid1 = 10,                # Axial points, tube 1
    Rgrid2 = 14,                # Radial points, tube 2
    Zgrid2 = 10,                # Axial points, tube 2

    # --- Numerics ----------------------------------------------------------------------
    model_mode = 'flowtube2',   # Chemistry (stiff solver, parallel over grid cells) and transport alternate
    dt = 2e-3,                  # Transport time step (s). The chemistry is solved separately with a stiff
                                # solver, so dt only has to keep diffusion/advection stable:
                                # dt < dr^2 / (2 D) and dt < dx / u_max. Finer grids need smaller dt
                                # (R = 60 used dt = 5e-5). NaN in the output = dt too large.
    timesteps = 250,            # Transport steps per iteration; chemistry and transport alternate every
                                # dt * timesteps = 0.5 s (same interval as the publication runs)
    inter_acc = 1e-3,           # Converged when IDHDP at the outlet changes by < 0.1 % between iterations
                                # (publication runs: 1e-5)
    fix_timstep = 2,            # Minimum number of iterations before convergence is accepted
    num_plot = 5,               # Print progress / update the plots every 5 iterations

    # --- Diffusion coefficients (cm2 s-1); the others are estimated from their formula
    Diff_setname = ['OH', 'HO2'],
    Diff_set     = [0.215, 0.141],

    # --- Chemical mechanism (input_mechanism/Isoprene/Wennberg/) ------------------------
    folder_mechaism = 'Isoprene/Wennberg',
    sch_name = 'isoprene_reduced_plus_v5_SUNfinal.eqn',  # Wennberg reduced mechanism + H2O2 photolysis (SUN)
    tsv_file = 'SMILE_formula.csv',     # Species SMILES table (species XML is chemical_species_custom.xml)
    flag_mech = '0',                    # '1' = MCM export, '0' = other mechanism

    # --- Species treatment -------------------------------------------------------------
    const_comp = ['O2', 'H2O'],         # Held constant everywhere (CO is added below)
    Init_comp = ['ISOP1OH2OOH'],        # Fixed at the tube inlet (from the ISOP1OH2OOHconc column) ...
    Init_set  = 'on',                   # ... 'on' = keep it fixed at the inlet. OH is NOT listed: it is
                                        # produced in the illuminated section by the itProd reaction
    key_spe_for_plot = 'IDHDP',         # Species used for the convergence test and progress output
    plot_spec = ['ISOP1OH2OOH', 'OH', 'HO2', 'IDHDP', 'ICPDH', 'IDHPE', 'HAC', 'GLYC',
                 'IHPOO1', 'IHPOO2', 'IHPOO3'],

    final_output_method = 'mean',       # Outlet concentration: 'mean' = area average over the cross-section
                                        # (default), 'weighted' = flow-weighted (what leaves the tube per
                                        # unit time, as in the Matlab model); see README section 7
    # --- Input / output ----------------------------------------------------------------
    file_name = 'Isoprene_continuousOH.csv',  # Input CSV in Input_files/; results in Export_files/Isoprene_continuousOH/

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

# Add a constant 100 ppb CO background (any extra constant species can be added this way)
CO_conc = 100e-9 * modelparams.p / (1.3806488e-23 * modelparams.TEMP * 1e6)   # molecules cm-3
modelparams.COconc = np.full_like(modelparams.O2conc, CO_conc)
modelparams.const_comp = list(modelparams.const_comp) + ['CO']
modelparams.const_comp_conc = np.transpose([getattr(modelparams, i + 'conc') for i in modelparams.const_comp])

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
