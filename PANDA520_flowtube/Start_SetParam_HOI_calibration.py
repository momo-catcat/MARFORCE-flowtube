# Example 2: HOI calibration — POINT OH source, I2 + OH -> HOI + I
#
# OH and HO2 are produced at the inlet by H2O photolysis and held there (Init_set = "on");
# the measured I2 (column I2conc of the input CSV) is held constant along the tube.
# Same point-source set-up as Example 1: no Rl / Ll, Itx = lamp It product (photons cm-2).
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

    # Tube geometry. Q1 == Q2, so both sections are simulated as one tube of length L1 + L2
    R1 = 0.78, # Inner radius of the 1st tube (cm)
    L1 = 41,  # Length of the 1st tube (cm)
    R2 = 1.04, # Inner radius of the 2nd tube (cm)
    L2 = 58.5,  # Length of the 2nd tube (cm)
    Itx = 4.84e10,  # Lamp It product at Qx (photons cm-2); see 10.5194/amt-4-437-2011
    Qx = 20,  # Flow (slpm) at which Itx was determined
    sampleflow = 22.5, # Inlet flow of the CIMS (slpm)

    outflowLocation = 'before',  # Outflow tube location: 'before' or 'after' injecting air, water, and I2
    fullOrSimpleModel = 'full',  # 'simple': Gormley & Kennedy approximation, 'full': flow model (slower)

    O2ratio = 0.209,  # O2 fraction in synthetic air

    # Grid (axial points = Zgrid1 + Zgrid2 over L1 + L2)
    Zgrid1 = 10,
    Rgrid1 = 20,
    Zgrid2 = 10,
    Rgrid2 = 20,

    dt = 1e-3, # Time step (s); flowtube1 integrates chemistry explicitly, so dt must resolve the fastest reaction
    timesteps = 500,  # Time steps per iteration (iteration length = dt * timesteps = 0.5 s)
    inter_acc = 1e-3, # Convergence: relative change of key species between iterations
    fix_timstep = 2, # Minimum number of iterations
    model_mode = 'flowtube1',  # Point source: chemistry solved together with transport at every time step

    num_plot = 5, # Plot/print every 5 iterations
    # Species with user-defined diffusion coefficients (cm2 s-1)
    Diff_setname = ['OH', 'HO2'],
    Diff_set = [0.215, 0.141],

    # Chemical mechanism (non-MCM: species/SMILES come from chemical_species_custom.xml in the same folder)
    sch_name = 'HOI_cali_chem.txt',
    tsv_file = 'chemical_species_custom.xml',
    folder_mechaism = 'HOI',
    flag_mech = '0', # 1: MCM, 0: other mechanism

    # Constant concentration species (I2 from the I2conc column of the input CSV)
    const_comp = ['I2', 'O2', 'H2O'],

    # Point OH source: OH and HO2 fixed at the inlet
    Init_comp = ['OH', 'HO2'],
    Init_set  = 'on',

    # Key species for the convergence criterion
    key_spe_for_plot = 'HOI',

    # Species to be plotted
    plot_spec = ['OH', 'HOI', 'HO2', 'I', 'I2'],

    # Input CSV in Input_files/ (one row per stage); also the output folder name
    file_name = 'HOI_calibration.csv',

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
