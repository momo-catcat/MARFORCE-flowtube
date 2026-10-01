# =====================================================================================
# Example 5: TRANSPORT-ONLY test (model_mode = 'kinetic') against Gormley & Kennedy
# =====================================================================================
# H2SO4 (1e8 cm-3) enters a 1 m tube (R = 0.78 cm) and is only transported: laminar flow,
# radial/axial diffusion and loss to the wall, no chemistry. The fraction that leaves the tube
# is compared with the analytical Gormley & Kennedy (1949) penetration for laminar flow,
#     mu = pi * D * L / Q,  P = 0.819 exp(-3.657 mu) + 0.097 exp(-22.3 mu) + 0.032 exp(-57 mu)
#     (P = 1 - 2.56 mu^(2/3) + 1.2 mu + 0.177 mu^(4/3) for mu < 0.02),
# which is the flow-weighted outlet / inlet ratio -> final_output_method = 'weighted'.
# Use this mode to check the transport core (grid, dt) of a new set-up before adding chemistry.
#
# Input: Input_files/Kinetic_test.csv (2 stages: 22.5 and 10 slpm).
# Run from PANDA520_flowtube/:   python Start_SetParam_kinetic_test.py
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

    # --- Tube geometry (one tube, no Rl/Ll) ------------------------------------------
    R1 = 0.78,                  # Inner RADIUS of the tube (cm)
    L1 = 100,                   # Length of the tube (cm)

    # --- Lamp (only needed so that every stage counts as active; no chemistry is run) -
    Itx = 4.84e10,              # Lamp It product (photons cm-2)
    Qx = 20,                    # Flow (slpm) at which Itx was measured

    # --- Flows ---------------------------------------------------------------------------
    sampleflow = 22.5,          # CIMS inlet flow (slpm)
    outflowLocation = 'before', # Exhaust 'before' or 'after' the gas injection
    O2ratio = 0.209,            # O2 fraction of the synthetic air

    # --- Grid (one tube: Zgrid1 + Zgrid2 axial points over L1) -------------------------
    Rgrid1 = 80,                # Radial grid points across the diameter
    Zgrid1 = 20,                # Axial grid points, part 1
    Rgrid2 = 80,                # Must equal Rgrid1
    Zgrid2 = 20,                # Axial grid points, part 2

    # --- Numerics ----------------------------------------------------------------------
    model_mode = 'kinetic',     # Transport only: laminar advection + diffusion + wall loss, NO chemistry
    dt = 2.5e-4,                # Time step (s); only transport stability matters (no chemistry)
    timesteps = 4000,           # Time steps per iteration; one iteration = 1 s (~1 residence time at 10 slpm)
    inter_acc = 1e-5,           # Converged when H2SO4 at the outlet changes by < 1e-5 between iterations
    fix_timstep = 3,            # Minimum number of iterations
    num_plot = 1,               # Print progress / update the plots every iteration

    # --- Diffusion coefficient of the test species (cm2 s-1) ----------------------------
    Diff_setname = ['H2SO4'],
    Diff_set     = [0.088],

    # --- Mechanism: only used for the species list (chemistry is switched off) ----------
    folder_mechaism = 'SA',
    sch_name = 'mcm_export.fac',
    tsv_file = 'mcm_export_species.tsv',
    flag_mech = '1',

    # --- Species treatment -------------------------------------------------------------
    const_comp = ['O2', 'H2O'],         # Held constant (no effect without chemistry)
    Init_comp = [],                     # The inlet species is set after calculate_concentrations() below
    Init_set  = 'on',                   # Keep the inlet species fixed at the inlet
    key_spe_for_plot = 'H2SO4',         # Species used for the convergence test
    plot_spec = ['H2SO4'],

    # --- Input / output ----------------------------------------------------------------
    final_output_method = 'weighted',   # Flow-weighted outlet value = what Gormley & Kennedy describe
    file_name = 'Kinetic_test.csv',     # Input CSV in Input_files/; results in Export_files/Kinetic_test/

    # --- Box-model and internal settings (not used here; keep as is) -------------------
    dil_fac_now = 0,
    wall_loss_set = 0,
    tot_time = 10,
    save_step = 0.5,
    pars_skip = 0,
)

'''''''''
Calculate the concentrations of every stage from the flows in the input CSV
'''''''''
calculate_concentrations(modelparams)

# Test species at the inlet: H2SO4 = 1e8 cm-3 in every stage
H2SO4_in = 1e8
modelparams.Init_comp = ['H2SO4']
modelparams.Init_comp_conc = np.full((len(modelparams.OHconc), 1), H2SO4_in)


def gormley_kennedy(D, L, Q):
    """Penetration through a tube in laminar flow (D cm2 s-1, L cm, Q cm3 s-1)."""
    mu = np.pi * D * L / Q
    if mu < 0.02:
        return 1 - 2.56 * mu ** (2 / 3) + 1.2 * mu + 0.177 * mu ** (4 / 3)
    return 0.819 * np.exp(-3.657 * mu) + 0.097 * np.exp(-22.3 * mu) + 0.032 * np.exp(-57 * mu)


'''''''''
Run the model and compare with Gormley & Kennedy
'''''''''
if __name__ == "__main__":
    import glob, os
    import pandas as pd
    multiprocessing.set_start_method("spawn")
    num_stage = modelparams.OHconc
    start_time = time.perf_counter()

    Run_flowtube(modelparams, num_stage)

    total_time = (time.perf_counter() - start_time) / 60
    if multiprocessing.current_process().name == "MainProcess":
        print(f"Total execution time: {total_time:.2f} mins")
        res_file = max(glob.glob(modelparams.export_file_folder + '*allstage_final_results.csv'), key=os.path.getmtime)
        res = pd.read_csv(res_file, index_col=0)
        print('\nTransport test: H2SO4 penetration (outlet / inlet)')
        print('  Q (slpm)    model    Gormley-Kennedy    difference')
        for j in range(len(res)):
            Q = float(modelparams.Q1[j]) / 60          # sccm -> cm3 s-1
            p_model = res['H2SO4'].iloc[j] / H2SO4_in
            p_gk = gormley_kennedy(0.088, float(modelparams.L1), Q)
            print(f'  {Q * 60 / 1000:8.1f}   {p_model:6.3f}        {p_gk:6.3f}        {100 * (p_model / p_gk - 1):+5.1f} %')
