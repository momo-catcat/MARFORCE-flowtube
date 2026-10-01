# Example 1: SA (H2SO4) calibration — POINT OH source, one tube
#
# OH and HO2 are produced at the tube inlet by H2O photolysis (UV lamp at the inlet):
#   [OH]0 = [HO2]0 = Itx * Qx / Q * 7.22e-20 * [H2O]
# They are held at the inlet (Init_set = "on") and react with SO2 down the tube to form H2SO4.
# To use a POINT source: do NOT set Rl / Ll, give Itx as the lamp It product (photons cm-2),
# list OH and HO2 in Init_comp and set Init_set = "on".
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
    TEMP = 298,  # Temperature, K

    # Tube geometry (one tube: no R2/L2). No Rl/Ll -> point OH source at the inlet
    R1 = 0.78, # Inner radius of the tube (cm)
    L1 = 26,  # Length of the tube (cm)
    Itx = 4.84e10,  # Lamp It product at Qx (photons cm-2); see 10.5194/amt-4-437-2011
    Qx = 20,  # Flow (slpm) at which Itx was determined
    sampleflow = 22.5, # Inlet flow of the CIMS (slpm), same as the total flow Q here

    outflowLocation = 'before',  # Outflow tube location: 'before' or 'after' injecting air, water, and SO2
    fullOrSimpleModel = 'full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (slower)

    SO2ratio = 5e-3,  # SO2 mole fraction in the gas bottle (5000 ppm)
    O2ratio = 0.209,  # O2 fraction in synthetic air

    # Grid (one tube: axial points = Zgrid1 + Zgrid2 over L1)
    Zgrid1 = 10,  # Axial grid points
    Rgrid1 = 20,  # Radial grid points
    Zgrid2 = 10,
    Rgrid2 = 20,

    dt = 1e-4, # Time step (s); flowtube1 integrates chemistry explicitly, so dt must resolve the fastest reaction
    timesteps = 2000,  # Time steps per iteration (iteration length = dt * timesteps = 0.2 s)
    inter_acc = 1e-3, # Convergence: relative change of key species between iterations
    fix_timstep = 2, # Minimum number of iterations
    model_mode = 'flowtube1',  # Point source: chemistry solved together with transport at every time step

    num_plot = 5, # Plot/print every 5 iterations
    # Species with user-defined diffusion coefficients (cm2 s-1)
    Diff_setname = ['OH', 'HO2', 'SO3', 'H2SO4', 'H', 'O'],
    Diff_set = [0.215, 0.141, 0.126, 0.088, 0.20, 0.20],

    # Chemical mechanism (MCM export + SMILES table)
    sch_name = 'mcm_export.fac',
    tsv_file = 'mcm_export_species.tsv',
    folder_mechaism = 'SA',
    flag_mech = '1', # 1: MCM, 0: other mechanism

    # Constant concentration species
    const_comp = ['SO2', 'O2', 'H2O'],

    # Point OH source: OH and HO2 fixed at the inlet
    Init_comp = ['OH', 'HO2'],
    Init_set  = 'on',

    # Key species for the convergence criterion
    key_spe_for_plot = 'H2SO4',

    # Species to be plotted
    plot_spec = ['OH', 'HO2', 'HSO3', 'SO3', 'H2SO4', 'H2O', 'SO2'],

    # Input CSV in Input_files/ (one row per stage); also the output folder name
    file_name = 'SA_calibration.csv',

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
