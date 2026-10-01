# MARFORCE-Flowtube

MARFORCE-Flowtube is a 2-D (radial × axial) solver for **convection–diffusion–reaction** processes in laminar flow reactors. It is used to simulate CIMS calibration set-ups (e.g. H2SO4 / OH production from H2O photolysis) and flow-tube oxidation experiments (e.g. isoprene + OH chemistry with the Wennberg mechanism), including:

- one tube, or two tubes of different diameter in series (with or without an extra Y-piece flow),
- a point OH source at the tube inlet **or** a continuous OH source along an illuminated section (UV lamp),
- arbitrary gas-phase mechanisms in MCM (`.fac`) or KPP (`.eqn`) format,
- a 0-D **box-model** mode for comparison with the flow-tube results,
- parallel (Numba + numbalsoda) chemistry integration, checkpointing and warm starts for long HPC runs.

The model was published with our manuscript in [Atmospheric Measurement Techniques](https://amt.copernicus.org/articles/16/4461/2023/) — please cite it if you use MARFORCE-Flowtube. Parts of the chemistry parser are adapted from [PyCHAM](https://github.com/simonom/PyCHAM) (Simon O'Meara, GPL-3.0).

Feedback is very welcome — please open an issue.

---

## Table of contents

1. [Repository layout](#1-repository-layout)
2. [Installation](#2-installation)
3. [Quick start — run an example](#3-quick-start)
4. [Setting up your own simulation](#4-setting-up-your-own-simulation)
   - [4.1 Input CSV (flows and measured concentrations)](#41-input-csv)
   - [4.2 Chemical mechanism](#42-chemical-mechanism)
   - [4.3 Start script parameters](#43-start-script-parameters)
   - [4.4 Model modes](#44-model-modes)
   - [4.5 Photolysis scaling (SUN)](#45-photolysis-scaling-sun)
5. [Outputs](#5-outputs)
6. [Running on an HPC cluster](#6-running-on-an-hpc-cluster)
7. [Tips and troubleshooting](#7-tips-and-troubleshooting)
8. [Acknowledgements](#8-acknowledgements)

---

## 1. Repository layout

```
MARFORCE-flowtube/
├── environment.yml / requirements.txt     # dependencies
└── PANDA520_flowtube/                     # run everything from this folder
    ├── Start_SetParam_SA.py               # example: H2SO4 / OH production (MCM, two tubes)
    ├── Start_SetParam_Isoprene.py         # example: isoprene + OH, Wennberg mechanism
    ├── Start_SetParam_Isoprene_SUN.py     # example: as above + H2O2/SUN photolysis (final mechanism)
    ├── Start_SetParam_box_model_comparison.py  # example: 0-D box model
    ├── run_flowtube_slurm.sh              # SLURM job template
    ├── Input_files/                       # input CSVs (one row = one experiment stage)
    ├── input_mechanism/                   # chemical mechanisms
    │   ├── SA/                            #   MCM SO2 → H2SO4 subset (.fac + species .tsv)
    │   ├── Isoprene/Wennberg/             #   Wennberg isoprene mechanism (.eqn, KPP format)
    │   ├── Isoprene/MCM/                  #   MCM isoprene subset
    │   └── HOI/                           #   HOI calibration mechanism (no example script yet)
    ├── photofiles/                        # MCM v3.2 cross-sections / quantum yields, lamp spectrum
    ├── Funcs/                             # model source code
    ├── kinetics/                          # diffusion coefficients and helper kinetics
    └── Export_files/                      # results are written here (one sub-folder per run)
```

> **Note:** `Funcs/rate_coeffs.py`, `Funcs/dydt_rec.py`, `Funcs/hyst_eq.py` and `Funcs/ode_solv.py` are **generated automatically** from the mechanism every time you run the model. Do not edit them by hand.

## 2. Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/momo-catcat/MARFORCE-flowtube.git
   cd MARFORCE-flowtube
   ```

2. Create the Python environment. Conda is recommended because **Open Babel** is easiest to install from conda-forge:

   ```bash
   conda env create -f environment.yml
   conda activate marforce
   ```

   Alternatively, with an existing Python ≥ 3.9:

   ```bash
   conda install -c conda-forge openbabel   # or your cluster's openbabel module
   pip install -r requirements.txt
   ```

   Required packages: `numpy`, `scipy`, `pandas`, `numba`, `numbalsoda`, `matplotlib`, `molmass`, `xmltodict`, `openbabel`. The code has been tested with Python 3.9/3.10, numpy 2.0, scipy 1.13, pandas 2.2, numba 0.60 and numbalsoda 0.3.4.

3. Check the installation:

   ```bash
   python -c "import numba, numbalsoda, molmass, xmltodict; from openbabel import pybel; print('OK')"
   ```

## 3. Quick start

Always run from inside `PANDA520_flowtube/` (the scripts import `Funcs.*` relative to it):

```bash
cd PANDA520_flowtube
python Start_SetParam_box_model_comparison.py     # ~30 s, good first test
python Start_SetParam_SA.py                       # H2SO4 / OH production in two tubes
python Start_SetParam_Isoprene.py                 # isoprene + OH, 13 experiment stages
python Start_SetParam_Isoprene_SUN.py             # final SUN-photolysis set-up (fine grid; use a cluster)
```

| Example script | Input CSV | Mechanism | Mode | Notes |
|---|---|---|---|---|
| `Start_SetParam_box_model_comparison.py` | `ISOP_OH_box_model_comparison.csv` | `Isoprene/Wennberg/isoprene_reduced_plus_v5_boxmodel.eqn` | `box` | No transport; OH, CO, ISOP1OH2OOH held constant. Runs in seconds. |
| `Start_SetParam_SA.py` | `OH_it_production.csv` | `SA/mcm_export.fac` (MCM) | `flowtube2` | OH from H2O photolysis in a lamp section → H2SO4. Two tubes. |
| `Start_SetParam_Isoprene.py` | `ISOP_OH_diff_OH.csv` | `Isoprene/Wennberg/isoprene_reduced_plus_v5.eqn` | `flowtube2` | Coarse grid (R=14). Several experiment stages. |
| `Start_SetParam_Isoprene_SUN.py` | `ISOP_OH_SUN.csv` | `Isoprene/Wennberg/isoprene_reduced_plus_v5_SUNfinal.eqn` | `flowtube2` | Final set-up: R=60 grid, `sun_a` photolysis scaling, 100 ppb CO background. Hours to days on one node. |

While running, the model prints the convergence of the key species for each iteration and shows/saves surface plots. Results go to `Export_files/<file_name without .csv>/` (see [Outputs](#5-outputs)).

On a machine without a display (e.g. a cluster) set `export MPLBACKEND=Agg`.

## 4. Setting up your own simulation

A simulation needs three things:

1. an **input CSV** in `Input_files/` describing the flows (and optionally measured concentrations) for every experiment stage,
2. a **chemical mechanism** in `input_mechanism/<folder>/`,
3. a **start script** `Start_SetParam_*.py` holding all model parameters.

The easiest way is to copy the example closest to your experiment and edit it.

### 4.1 Input CSV

Each row is one **experiment stage** (e.g. a different H2O flow or lamp setting); the model simulates each stage in turn. Column names are parsed by suffix:

- `<X>flow` — a flow in **sccm**. A trailing digit selects the tube: `H2Oflow1` = first tube, `H2Oflow2` = flow added at the Y-piece / second tube. Without a digit it is a single-tube set-up.
- `<X>conc` — a measured or prescribed concentration in **molecules cm⁻³** (e.g. `ISOP1OH2OOHconc`, `COconc`, `H2Oconc`). It becomes available as `modelparams.<X>conc` and can be used in `const_comp` or `Init_comp`.
- `Q` (single tube) or `Q1` and `Q2` (two tubes) — total flow in each tube, **sccm**. `Q2` is the flow in the second tube, *not* the Y-piece flow.
- `T` *(optional)* — temperature per stage (K); otherwise `TEMP` from the start script is used.
- `time` *(optional)* — residence/run time per stage (s), used to extend the minimum number of iterations.
- Other columns (e.g. `UVC`, `Start time`) are ignored.

The tube configuration is detected automatically:

| Columns present | Set-up |
|---|---|
| `H2Oflow`, `N2flow`, `Q`, no `L2` in start script | one tube |
| `H2Oflow`, `N2flow`, `Q` and `R2`, `L2` in start script | two tubes of different diameter, same flow |
| `H2Oflow1`, `N2flow1`, `Q1`, `Q2`, … | two tubes, extra flow added at the Y-piece |

**Gas-bottle species.** For every `<X>ratio` parameter in the start script (the mole fraction of X in its gas bottle, e.g. `SO2ratio = 5e-3`) and a matching `<X>flow` column, the concentration is computed from the flow dilution. `O2` is computed from `O2flow` × `O2ratio`, and `H2O` from `H2Oflow` assuming saturated water vapour (measured `H2Oconc` overrides this when the H2O flow is ≥ 1000 sccm).

**OH production.** `[OH]₀ = (Itx · Qx / Q) · σ(H2O, 185 nm) · [H2O]`, with σ = 7.22×10⁻²⁰ cm². `Itx` is the lamp's *It* product measured at flow `Qx` (see [doi:10.5194/amt-4-437-2011](https://doi.org/10.5194/amt-4-437-2011)). If `Rl > 0` in the start script the OH source is treated as **continuous** along the illuminated section (`Rl`, `Ll`); otherwise OH is a **point** source at the inlet.

Examples: `Input_files/OH_it_production.csv` (two tubes + SO2), `Input_files/SA_cali_2021-09-10.csv` (one tube), `Input_files/ISOP_OH_diff_OH.csv` (two tubes + prescribed ISOP1OH2OOH, 13 stages).

### 4.2 Chemical mechanism

Put the mechanism in its own folder under `input_mechanism/` and point the start script to it with `folder_mechaism` and `sch_name`.

| Format | `flag_mech` | `sch_name` | `tsv_file` |
|---|---|---|---|
| MCM export (`.fac`, KPP-style from the [MCM website](https://mcm.york.ac.uk/MCM)) | `'1'` | e.g. `mcm_export.fac` | the MCM species export `mcm_export_species.tsv` (SMILES); `chemical_species_custom.xml` is generated from it automatically |
| Other KPP-style mechanisms (`.eqn`, e.g. Wennberg isoprene) | `'0'` | e.g. `isoprene_reduced_plus_v5.eqn` | a CSV of species → SMILES (`SMILE_formula.csv`); you must also provide `chemical_species_custom.xml` in the same folder |

Rate expressions may use Fortran or Python scientific notation and the functions `EXP/exp, dsqrt, dlog, LOG, dabs, LOG10, numpy.*`, and may depend on `TEMP` (K), `RH` (0–1), `M`, `N2`, `O2` (molecules cm⁻³). Additional helpers are available: Wennberg-style `TROE, TUN, ALK, NIT, EPO, ISO1, ISO2, KCO`, `itProd` (lamp *It* product), `SUN` (photolysis scaling, see [4.5](#45-photolysis-scaling-sun)) and the MCM photolysis rates `J(n)`.

### 4.3 Start script parameters

All parameters live in one `SimulationParams(...)` call at the top of each start script. Numbers are converted to `float32` automatically, grid sizes to `int`.

**Experimental conditions**

| Parameter | Description |
|---|---|
| `p` | Pressure (Pa) |
| `TEMP` | Temperature (K) |
| `R1`, `L1` | Inner radius (cm) and length (cm) of the first tube |
| `R2`, `L2` | Radius and length of the second tube (omit `L2` for a single tube) |
| `Rl`, `Ll` | Radius and length of the illuminated (lamp) section; `Rl > 0` ⇒ continuous OH source |
| `Itx`, `Qx` | Lamp *It* product and the flow (slpm) at which it was determined |
| `sampleflow` | CIMS inlet flow (slpm) |
| `outflowLocation` | `'before'` or `'after'` — whether the exhaust is before or after the injection of air/H2O/reactant |
| `fullOrSimpleModel` | `'full'` flow model or `'simple'` Gormley–Kennedy approximation |
| `O2ratio`, `<X>ratio` | O2 fraction in synthetic air; bottle mixing ratio of species X |
| `sun_a` *(optional)* | Value of `SUN` in the mechanism (photolysis scaling); default 0 |

**Numerics**

| Parameter | Description |
|---|---|
| `Rgrid1`, `Zgrid1` / `Rgrid2`, `Zgrid2` | Radial / axial grid points for tube 1 / tube 2 |
| `Rgridl`, `Zgridl` | Grid of the illuminated section (continuous OH source) |
| `dt` | Transport time step (s). Reduce it if you see `NaN` or oscillations (finer grids need smaller `dt`). |
| `timesteps` | Transport steps per iteration |
| `inter_acc` | Convergence threshold on the relative change of `key_spe_for_plot` between iterations |
| `fix_timstep` | Minimum number of iterations before convergence is checked |
| `num_plot` | Plot/print every *n* iterations |
| `model_mode` | See [4.4](#44-model-modes) |
| `Diff_setname`, `Diff_set` | Species with user-defined diffusion coefficients (cm² s⁻¹); others are computed from their formula |
| `tot_time`, `save_step` | Box-model integration time and output interval (s) |
| `dil_fac_now`, `wall_loss_set` | Box-model dilution rate (s⁻¹) and wall loss (0 = off, 1 = computed, other = constant value) |
| `pars_skip` | 0 = parse mechanism (normal); 1 = skip parsing |

**Chemistry and outputs**

| Parameter | Description |
|---|---|
| `folder_mechaism`, `sch_name`, `tsv_file`, `flag_mech` | Mechanism location and format (see [4.2](#42-chemical-mechanism)) |
| `const_comp` | Species held constant everywhere in the tube (e.g. `['SO2','O2','H2O']`). Their concentrations must exist as `modelparams.<X>conc` (from the CSV or computed). |
| `Init_comp`, `Init_set` | Species set at the first axial grid point (inlet); `Init_set='on'` keeps them fixed there |
| `key_spe_for_plot` | Species used for the convergence criterion and progress output (e.g. `'H2SO4'`, `'IDHDP'`) |
| `plot_spec` | Species to plot |
| `file_name` | Input CSV name in `Input_files/`; also the name of the output sub-folder |
| `input_file_folder`, `export_file_folder`, `input_mechanism_folder` *(optional)* | Override the default folders |

After the parameter block, `calculate_concentrations(modelparams)` turns the CSV into concentrations. You can modify concentrations afterwards — e.g. `Start_SetParam_Isoprene_SUN.py` adds a constant 100 ppb CO background:

```python
calculate_concentrations(modelparams)
CO_conc = 100e-9 * modelparams.p / (1.3806488e-23 * modelparams.TEMP * 1e6)
modelparams.COconc = np.full_like(modelparams.O2conc, CO_conc)
modelparams.const_comp = list(modelparams.const_comp) + ['CO']
modelparams.const_comp_conc = np.transpose([getattr(modelparams, i + 'conc') for i in modelparams.const_comp])
```

The number of stages to run is `num_stage = modelparams.OHconc` (all rows with OH > 0); set it to an integer *N* to run only the first *N* stages.

### 4.4 Model modes

| `model_mode` | Description |
|---|---|
| `'flowtube1'` | Flow tube with small mechanisms (≲ 50 reactions), point OH source |
| `'flowtube2'` | Flow tube with large mechanisms and/or a continuous OH source (parallel Numba + numbalsoda solver) — used in all flow-tube examples |
| `'box'` | 0-D box model: chemistry only, no diffusion or advection |
| `'kinetic'` | Transport only, no chemistry (tests the convection–diffusion core) |

### 4.5 Photolysis scaling (SUN)

Mechanisms can contain photolysis reactions written as `SUN * k`, e.g. in `isoprene_reduced_plus_v5_SUNfinal.eqn`:

```
% SUN*5e-5: H2O2 = OH + OH;
```

`SUN` takes the value of `sun_a` from the start script (`sun_a = 0` switches these reactions off). `Start_SetParam_Isoprene_SUN.py` uses `sun_a = 2e5`; to scan photolysis strength, copy the script and change `sun_a` and `file_name` (and the CSV name) for each case.

## 5. Outputs

Results are written to `Export_files/<file_name without .csv>/`. File names encode the run settings, e.g. `R60L30inter_acc1e-06dt5e-05timstep1e+04itx1.9e-10flowtube2_...`:

| File | Content |
|---|---|
| `*_allstage_final_results.csv` | Mean concentration (molecules cm⁻³) at the tube exit of **every species** (columns) for each stage (rows) — the main result |
| `*_Stage_<j>_endtube_profile_output.txt` | Radial profile of all species at the tube end for stage *j* |
| `*_Stage_<j>_1st_tube_profile_output.txt` | Same, at the end of the first tube (two-tube set-ups) |
| `*_params.csv` | All parameters used for the run (for reproducibility) |
| `*_Stage_<j>_delta_keyspecforplot.csv` | Convergence history of `key_spe_for_plot` (flow-tube modes); box mode writes `<file_name>R…box.csv` |
| `warmstart_stage<j>_*.npz` | Converged 2-D field, used to warm-start the next stage or a re-run |
| `*checkpoint*.npz` | Periodic checkpoint; a re-run with the same settings **resumes** from it automatically |

Delete the `*.npz` files in the output folder if you want a fresh start.

## 6. Running on an HPC cluster

Fine grids (e.g. `Rgrid = 60`) and large mechanisms are expensive. `run_flowtube_slurm.sh` is a SLURM template; edit the account and environment lines, then:

```bash
cd PANDA520_flowtube
mkdir -p logs
sbatch run_flowtube_slurm.sh Start_SetParam_Isoprene_SUN.py
```

The chemistry step is parallelised over grid cells with Numba — set `NUMBA_NUM_THREADS` to the number of CPUs (done in the template) and keep `OMP/MKL/OPENBLAS_NUM_THREADS=1`. For parameter scans, make one start script per case (different `file_name`) and submit them as a job array. If several jobs share the same folder, give each its own copy of the code (the auto-generated `Funcs/*.py` files are rewritten at the start of each run).

## 7. Tips and troubleshooting

- **`NaN` or exploding concentrations** — reduce `dt` (e.g. 1e-4 → 5e-5 for `Rgrid = 60`), or coarsen the grid.
- **No convergence / too slow** — start with a coarse grid (`Rgrid ≈ 14–30`) and a looser `inter_acc`, then refine. Checkpoints and warm starts let you continue interrupted runs.
- **`ModuleNotFoundError: Funcs`** — run the scripts from inside `PANDA520_flowtube/` (or add it to `PYTHONPATH`).
- **Open Babel import error** — install it with `conda install -c conda-forge openbabel` (pip wheels are often unavailable).
- **Species not found** — names in `const_comp`, `Init_comp`, `plot_spec`, `key_spe_for_plot` and `Diff_setname` must match the mechanism exactly; species in `const_comp`/`Init_comp` need a concentration (`<X>conc` column or `<X>ratio` + `<X>flow`).
- **Multiprocessing on macOS/Windows** — keep the `if __name__ == "__main__":` block and `multiprocessing.set_start_method("spawn")` from the examples.

## 8. Acknowledgements

We thank the ACCC Flagship funded by the Academy of Finland grant number 337549, Academy professorship funded by the Academy of Finland (grant no. 302958), Academy of Finland projects no. 346370, 325656, 316114, 314798, 325647, 341349 and 349659. European Research Council (ERC) project ATM-GTP Contract No. 742206. The Arena for the gap analysis of the existing Arctic Science Co-Operations (AASCO) funded by Prince Albert Foundation Contract No 2859. M.K. thanks the Jane and Aatos Erkko Foundation for providing funding. M.K. and X.-C.H thank the Jenny and Antti Wihuri Foundation for providing funding for this research.
