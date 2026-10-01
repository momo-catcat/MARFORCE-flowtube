# =====================================================================================
# Example 1b: H2SO4 (SA) calibration, POINT OH source, ODE-solver chemistry (flowtube2)
# =====================================================================================
# A UV lamp at the tube inlet photolyses H2O; the OH and HO2 formed there,
#     [OH]0 = [HO2]0 = Itx * Qx / Q * 7.22e-20 cm2 * [H2O],
# are held at the inlet (Init_set = 'on') and react with SO2 along the tube:
#     OH + SO2 -> HSO3,  HSO3 + O2 -> HO2 + SO3,  SO3 + 2 H2O -> H2SO4.
# The model returns the mean H2SO4 at the tube outlet, which is what the CIMS samples.
#
# Point-source set-up = no Rl/Ll, Itx is the lamp It product, Init_comp = ['OH','HO2'],
# Init_set = 'on'.
# Same case as Example 1a, but the chemistry is solved with the stiff ODE solver (LSODA, parallel over
# grid cells), alternating with transport every dt * timesteps = 1 ms. This interval must be much shorter
# than the residence time (~0.13 s); both examples agree within 0.2 % (README section 3).
# Two tube sections of different radius (R2/L2) or a Y-piece (two-flow CSV) can be added: see README 5.1.
# Input: Input_files/SA_calibration.csv (2 stages).
# Run from PANDA520_flowtube/:   python Start_SetParam_SA_calibration_ODE.py
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
    TEMP = 298,                 # Temperature (K); a 'T' column in the input CSV overrides it per stage

    # --- Tube geometry (one tube; no R2/L2). No Rl/Ll -> point OH source -------------
    R1 = 0.78,                  # Inner RADIUS of the tube (cm)
    L1 = 26,                    # Length of the tube from the lamp to the CIMS inlet (cm)

    # --- OH source (point) -------------------------------------------------------------
    Itx = 4.84e10,              # Lamp It product (photons cm-2) measured at flow Qx; see doi:10.5194/amt-4-437-2011
    Qx = 20,                    # Flow (slpm) at which Itx was measured; It is scaled to the actual flow Q

    # --- Flows and gas bottles ---------------------------------------------------------
    sampleflow = 22.5,          # CIMS inlet flow (slpm); equal to the tube flow Q in this one-tube set-up
    outflowLocation = 'before', # Exhaust 'before' or 'after' the gas injection: dilution uses Q ('before')
                                # or the sum of the injected flows ('after')
    SO2ratio = 5e-3,            # SO2 mole fraction in the gas bottle (5e-3 = 5000 ppm); used with the SO2flow column
    O2ratio = 0.209,            # O2 fraction of the synthetic air in the O2flow column

    # --- Grid (single flow: the axial grid is Zgrid1 + Zgrid2 points over L1) ---------
    Rgrid1 = 80,                # Radial grid points across the diameter (80 = original MARFORCE resolution)
    Zgrid1 = 20,                # Axial grid points, part 1
    Rgrid2 = 80,                # Must equal Rgrid1
    Zgrid2 = 20,                # Axial grid points, part 2 (40 axial points in total)

    # --- Numerics ----------------------------------------------------------------------
    model_mode = 'flowtube2',   # ODE: chemistry (stiff solver) and transport alternate
    dt = 2.5e-4,                # Transport time step (s). The chemistry is solved by the stiff ODE solver, so
                                # dt only has to keep diffusion/advection stable; 2.5e-4 s is stable on this
                                # R80 grid (reduce it if NaN appears, e.g. on a finer grid)
    timesteps = 4,              # Transport steps per iteration: chemistry and transport alternate every
                                # dt * timesteps = 1 ms (<< residence time 0.13 s; 0.5 or 2 ms give the same result)
    inter_acc = 1e-6,           # Converged when H2SO4 at the outlet changes by < 1e-6 between 1-ms iterations
    fix_timstep = 400,          # Minimum number of iterations (400 x 1 ms = 0.4 s, ~3 residence times)
    num_plot = 100,             # Print progress / update the plots every 100 iterations

    # --- Diffusion coefficients (cm2 s-1) for selected species; the others are estimated from their formula
    Diff_setname = ['OH', 'HO2', 'SO3', 'H2SO4', 'H', 'O'],
    Diff_set     = [0.215, 0.141, 0.126, 0.088, 0.20, 0.20],

    # --- Chemical mechanism (input_mechanism/SA/) ---------------------------------------
    folder_mechaism = 'SA',                 # Sub-folder of input_mechanism/
    sch_name = 'mcm_export.fac',            # Mechanism file (MCM export format)
    tsv_file = 'mcm_export_species.tsv',    # MCM species list with SMILES (used to build the species XML)
    flag_mech = '1',                        # '1' = MCM export, '0' = other mechanism

    # --- Species treatment -------------------------------------------------------------
    const_comp = ['SO2', 'O2', 'H2O'],  # Held constant everywhere (excess reactants)
    Init_comp = ['OH', 'HO2'],          # Fixed at the tube inlet (point OH source) ...
    Init_set  = 'on',                   # ... 'on' = keep them fixed at the inlet
    key_spe_for_plot = 'H2SO4',         # Species used for the convergence test and progress output
    plot_spec = ['OH', 'HO2', 'HSO3', 'SO3', 'H2SO4', 'H2O', 'SO2'],  # Species shown in the plots

    final_output_method = 'mean',       # Outlet concentration: 'mean' = area average over the cross-section
                                        # (default), 'weighted' = flow-weighted (what leaves the tube per
                                        # unit time, as in the Matlab model); see README section 7
    # --- Input / output ----------------------------------------------------------------
    file_name = 'SA_calibration.csv',   # Input CSV in Input_files/ (one row = one stage);
                                        # results go to Export_files/SA_calibration/

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
