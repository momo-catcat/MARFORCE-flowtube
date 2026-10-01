# MARFORCE-Flowtube

MARFORCE-Flowtube is a 2-D (radial × axial) solver for **convection–diffusion–reaction** processes in laminar flow reactors. It is used to simulate CIMS calibration set-ups (H2SO4 and HOI calibration with OH from H2O photolysis) and flow-tube oxidation experiments (e.g. isoprene chemistry with the Wennberg mechanism), including:

- one tube, or two tubes of different diameter in series (with or without an extra Y-piece flow),
- a **point** OH source at the tube inlet **or** a **continuous** OH source along an illuminated section (UV lamp),
- arbitrary gas-phase mechanisms in MCM (`.fac`) or KPP (`.eqn`) format, including photolysis scaling (`SUN`),
- a 0-D **box-model** mode,
- parallel (Numba + numbalsoda) chemistry integration, checkpointing and warm starts for long HPC runs.

The model was published with our manuscript in [Atmospheric Measurement Techniques](https://amt.copernicus.org/articles/16/4461/2023/) — please cite it if you use MARFORCE-Flowtube. Parts of the chemistry parser are adapted from [PyCHAM](https://github.com/simonom/PyCHAM) (Simon O'Meara, GPL-3.0).

Feedback is very welcome — please open an issue.

---

## Table of contents

1. [Repository layout](#1-repository-layout)
2. [Installation](#2-installation)
3. [Quick start — the four examples](#3-quick-start)
4. [Setting up your own simulation](#4-setting-up-your-own-simulation)
   - [4.1 Input CSV](#41-input-csv)
   - [4.2 OH source: point or continuous](#42-oh-source-point-or-continuous)
   - [4.3 Chemical mechanism](#43-chemical-mechanism)
   - [4.4 Start script parameters](#44-start-script-parameters)
   - [4.5 Model modes](#45-model-modes)
   - [4.6 Photolysis scaling (SUN)](#46-photolysis-scaling-sun)
   - [4.7 Speed vs. accuracy](#47-speed-vs-accuracy)
5. [Outputs](#5-outputs)
6. [Running on an HPC cluster](#6-running-on-an-hpc-cluster)
7. [Tips and troubleshooting](#7-tips-and-troubleshooting)
8. [Acknowledgements](#8-acknowledgements)
9. [License](#9-license)

---

## 1. Repository layout

```
MARFORCE-flowtube/
├── LICENSE                                     # GNU GPL v3
├── environment.yml / requirements.txt          # dependencies
└── PANDA520_flowtube/                          # run everything from this folder
    ├── Start_SetParam_SA_calibration.py        # Example 1: SA calibration, point OH source
    ├── Start_SetParam_HOI_calibration.py       # Example 2: HOI calibration, point OH source
    ├── Start_SetParam_Isoprene_box.py          # Example 3: isoprene chemistry, box model
    ├── Start_SetParam_Isoprene_continuousOH.py # Example 4: isoprene, continuous OH source + SUN photolysis
    ├── run_flowtube_slurm.sh                   # SLURM job template
    ├── Input_files/                            # input CSVs (one row = one experiment stage)
    ├── input_mechanism/                        # chemical mechanisms
    │   ├── SA/                                 #   MCM SO2 → H2SO4 subset (.fac + species .tsv)
    │   ├── HOI/                                #   iodine / HOI calibration mechanism
    │   ├── Isoprene/Wennberg/                  #   Wennberg isoprene mechanism (.eqn, KPP format)
    │   └── Isoprene/MCM/                       #   MCM isoprene subset
    ├── photofiles/                             # MCM v3.2 cross-sections / quantum yields, lamp spectrum
    ├── Funcs/                                  # model source code
    ├── kinetics/                               # diffusion coefficients and helper kinetics
    └── Export_files/                           # results are written here (one sub-folder per run)
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

   Required packages: `numpy`, `scipy`, `pandas`, `numba`, `numbalsoda`, `matplotlib`, `molmass`, `xmltodict`, `openbabel`. Tested with Python 3.9/3.10, numpy 2.0, scipy 1.13, pandas 2.2, numba 0.60 and numbalsoda 0.3.4.

3. Check the installation:

   ```bash
   python -c "import numba, numbalsoda, molmass, xmltodict; from openbabel import pybel; print('OK')"
   ```

## 3. Quick start

Always run from inside `PANDA520_flowtube/` (the scripts import `Funcs.*` relative to it):

```bash
cd PANDA520_flowtube
python Start_SetParam_SA_calibration.py          # Example 1, ~0.5 min
python Start_SetParam_HOI_calibration.py         # Example 2, ~0.5 min
python Start_SetParam_Isoprene_box.py            # Example 3, ~0.5 min
python Start_SetParam_Isoprene_continuousOH.py   # Example 4, ~4 min
```

The run times were measured on a laptop (10 cores). Each flow-tube example simulates **two experiment stages** (two rows of its input CSV) on a coarse grid, so that it runs quickly. For publication-quality results refine the grid (see [4.7](#47-speed-vs-accuracy)).

| # | Script | Input CSV | Mechanism | OH source | Mode | What it shows |
|---|---|---|---|---|---|---|
| 1 | `Start_SetParam_SA_calibration.py` | `SA_calibration.csv` | `SA/mcm_export.fac` (MCM) | point | `flowtube1` | Classic H2SO4 calibration: OH/HO2 from H2O photolysis at the inlet react with SO2 → H2SO4. One tube, 2 H2O-flow stages. |
| 2 | `Start_SetParam_HOI_calibration.py` | `HOI_calibration.csv` | `HOI/HOI_cali_chem.txt` | point | `flowtube1` | HOI calibration: inlet OH reacts with measured I2 (`I2conc` column, held constant) → HOI. 2 stages. |
| 3 | `Start_SetParam_Isoprene_box.py` | `Isoprene_box.csv` | `Isoprene/Wennberg/isoprene_reduced_plus_v5_SUNfinal.eqn` | held constant | `box` | 0-D chamber-style run with dilution and species-specific wall losses. Box mode uses one stage (first CSV row). |
| 4 | `Start_SetParam_Isoprene_continuousOH.py` | `Isoprene_continuousOH.csv` | same as Example 3 | continuous | `flowtube2` | ISOP1OH2OOH + OH in two tubes, OH produced along the lamp section, H2O2 photolysis via `sun_a` (SUN), 100 ppb CO background. 2 stages. |

While running, the model prints the convergence of the key species and shows surface plots. On a machine without a display set `export MPLBACKEND=Agg`. Results go to `Export_files/<input CSV name>/` (see [Outputs](#5-outputs)).

## 4. Setting up your own simulation

A simulation needs three things:

1. an **input CSV** in `Input_files/` with the flows (and optionally measured concentrations) of every experiment stage,
2. a **chemical mechanism** in `input_mechanism/<folder>/`,
3. a **start script** `Start_SetParam_*.py` holding all model parameters.

The easiest way is to copy the example closest to your experiment, rename the input CSV and edit the parameters.

### 4.1 Input CSV

Each row is one **experiment stage** (e.g. a different H2O flow); the model simulates the stages one after another. Only include the stages you want to run. Column names are parsed by suffix:

- `<X>flow` — a flow in **sccm**. A trailing digit selects the tube: `H2Oflow1` = first tube, `H2Oflow2` = flow added at the Y-piece / second tube. Without a digit it is a single-flow set-up.
- `<X>conc` — a measured or prescribed concentration in **molecules cm⁻³** (e.g. `I2conc`, `ISOP1OH2OOHconc`, `H2Oconc`). It becomes `modelparams.<X>conc` and can be used in `const_comp` or `Init_comp`.
- `Q` (single flow) or `Q1` and `Q2` (two flows) — total flow in each tube, **sccm**. `Q2` is the flow in the second tube, *not* the Y-piece flow.
- `T` *(optional)* — temperature per stage (K); otherwise `TEMP` from the start script is used.
- `time` *(optional)* — residence/run time per stage (s), used to extend the minimum number of iterations.
- Other columns (e.g. `UVC`, timestamps) are ignored.

The tube configuration is detected automatically:

| Columns / parameters | Set-up | Example |
|---|---|---|
| `H2Oflow`, `N2flow`, `Q`; no `L2` in the start script | one tube | 1 |
| `H2Oflow`, `N2flow`, `Q`; `R2`, `L2` in the start script | two tube sections, same flow (simulated as one tube of length `L1 + L2`) | 2 |
| `H2Oflow1`, `N2flow1`, `Q1`, `Q2`, … | two tubes, extra flow added at the Y-piece | 4 |

**Gas-bottle species.** For every `<X>ratio` parameter in the start script (mole fraction of X in its gas bottle, e.g. `SO2ratio = 5e-3`) and a matching `<X>flow` column, the concentration is computed from the flow dilution. `O2` is computed from `O2flow` × `O2ratio`, and `H2O` from `H2Oflow` assuming saturated water vapour (a measured `H2Oconc` column overrides this when the H2O flow is ≥ 1000 sccm).

### 4.2 OH source: point or continuous

This is the most important choice of a set-up. The model decides it from `Rl`: **if `Rl > 0` the source is continuous, otherwise it is a point source.**

| | **Point source** (Examples 1, 2) | **Continuous source** (Example 4) |
|---|---|---|
| Physical picture | OH/HO2 are made in a short lamp region at the inlet and then react downstream (classic SA / HOI calibration) | OH is produced all along an illuminated section of the tube while the chemistry proceeds |
| Geometry parameters | `R1`, `L1` (+ `R2`, `L2`); **do not set `Rl`, `Ll`** | `Rl`, `Ll`, `Rgridl`, `Zgridl` for the illuminated section, plus `R1`, `L1`, `R2`, `L2` |
| `Itx` | lamp *It* product at `Qx` (photons cm⁻², e.g. `4.84e10`) | first-order rate coefficient of `H2O = OH + H` in the lamp section (s⁻¹, e.g. `1.9e-10`) |
| How OH enters | `[OH]₀ = [HO2]₀ = Itx · Qx / Q · 7.22×10⁻²⁰ · [H2O]` is fixed at the inlet ([H2O] = measured `H2Oconc` if given and H2O flow ≥ 1000 sccm, otherwise calculated from the flow) | through the mechanism reaction `% itProd(modelparams) : H2O = OH + H ;` (must be in the mechanism) |
| `Init_comp` / `Init_set` | `Init_comp = ['OH', 'HO2']`, `Init_set = 'on'` | OH not in `Init_comp`; use `Init_set = 'on'` only for other species fixed at the inlet (e.g. `['ISOP1OH2OOH']`), otherwise `'off'` |
| `model_mode` | `'flowtube1'` | `'flowtube2'` |

To turn a continuous-source script into a point-source calibration: delete `Rl`, `Ll`, `Rgridl`, `Zgridl`; set `Itx` to the lamp *It* product; set `Init_comp = ['OH', 'HO2']`, `Init_set = 'on'`, `model_mode = 'flowtube1'` and choose `dt` small enough for the chemistry (see [4.7](#47-speed-vs-accuracy)). Example 1 is exactly this set-up.

### 4.3 Chemical mechanism

Put the mechanism in its own folder under `input_mechanism/` and point the start script to it with `folder_mechaism` and `sch_name`.

| Format | `flag_mech` | `sch_name` | `tsv_file` |
|---|---|---|---|
| MCM export (`.fac`, from the [MCM website](https://mcm.york.ac.uk/MCM)) | `'1'` | e.g. `mcm_export.fac` | the MCM species export `mcm_export_species.tsv` (SMILES); `chemical_species_custom.xml` is generated from it automatically |
| Other mechanisms (`.eqn` / `.txt`, e.g. Wennberg isoprene, HOI) | `'0'` | e.g. `isoprene_reduced_plus_v5_SUNfinal.eqn` | any file name; you must provide `chemical_species_custom.xml` (species → SMILES) in the same folder |

Mechanism file format (see `input_mechanism/HOI/HOI_cali_chem.txt` for a short example):

```
* comment lines start with * ;
KMT06 = 1 + (1.40D-21*EXP(2200/TEMP)*H2O) ;          <- generic rate coefficient (name < 10 characters)
% 2.1D-10 : I2 + OH = HOI + I ; # reference           <- reaction: % rate : reactants = products ;
```

Rate expressions may use Fortran or Python scientific notation and the functions `EXP/exp, dsqrt, dlog, LOG, dabs, LOG10, numpy.*`, and may depend on `TEMP` (K), `p` (Pa), `RH` (0–1), `M`, `N2`, `O2`, `H2O` (molecules cm⁻³). Additional helpers: Wennberg-style `TROE, TUN, ALK, NIT, EPO, ISO1, ISO2, KCO`, `itProd(modelparams)` (continuous OH source), `SUN` (see [4.6](#46-photolysis-scaling-sun)), `kDimer/kTrimer` (H2SO4 clustering) and the MCM photolysis rates `J(n)`. Do not leave blank lines in the mechanism file.

### 4.4 Start script parameters

All parameters live in one `SimulationParams(...)` call at the top of each start script. Numbers are converted to `float32` automatically, grid sizes to `int`.

**Experimental conditions**

| Parameter | Description |
|---|---|
| `p`, `TEMP` | Pressure (Pa) and temperature (K) |
| `R1`, `L1` | Inner radius (cm) and length (cm) of the first tube |
| `R2`, `L2` | Radius and length of the second tube (omit `L2` for a single tube) |
| `Rl`, `Ll` | Radius and length of the illuminated section — only for a continuous OH source ([4.2](#42-oh-source-point-or-continuous)) |
| `Itx`, `Qx` | Lamp parameter and the flow (slpm) at which it was determined ([4.2](#42-oh-source-point-or-continuous)) |
| `sampleflow` | CIMS inlet flow (slpm) |
| `outflowLocation` | `'before'` or `'after'` — whether the exhaust is before or after the injection of air/H2O/reactant |
| `fullOrSimpleModel` | `'full'` flow model or `'simple'` Gormley–Kennedy approximation |
| `O2ratio`, `<X>ratio` | O2 fraction in synthetic air; bottle mole fraction of species X |
| `sun_a` *(optional)* | Value of `SUN` in the mechanism ([4.6](#46-photolysis-scaling-sun)); default 0 |

**Numerics**

| Parameter | Description |
|---|---|
| `Rgrid1`, `Zgrid1` / `Rgrid2`, `Zgrid2` | Radial / axial grid points for tube 1 / tube 2 (single flow: axial points = `Zgrid1 + Zgrid2` over `L1 + L2`) |
| `Rgridl`, `Zgridl` | Grid of the illuminated section (continuous source only) |
| `dt` | Time step (s). Reduce it if you see `NaN` |
| `timesteps` | Time steps per iteration; one iteration covers `dt × timesteps` seconds |
| `inter_acc` | Convergence threshold: relative change of `key_spe_for_plot` at the tube end between two iterations |
| `fix_timstep` | Minimum number of iterations before convergence is accepted |
| `num_plot` | Plot/print every *n* iterations |
| `model_mode` | See [4.5](#45-model-modes) |
| `Diff_setname`, `Diff_set` | Species with user-defined diffusion coefficients (cm² s⁻¹); others are computed from their formula |
| `tot_time`, `save_step` | Box-model integration time and output interval (s) |
| `dil_fac_now`, `wall_loss_set`, `wall_loss_custom` | Box-model dilution rate (s⁻¹), wall loss (0 = off, 1 = computed, other = constant value) and per-species wall loss (s⁻¹) |
| `pars_skip` | 0 = parse mechanism (normal); 1 = skip parsing |
| `use_restart` *(optional)* | Checkpoints, warm start and stage carry-over ([5](#5-outputs)). Default: on for a continuous OH source, off for a point source |

**Chemistry and outputs**

| Parameter | Description |
|---|---|
| `folder_mechaism`, `sch_name`, `tsv_file`, `flag_mech` | Mechanism location and format ([4.3](#43-chemical-mechanism)) |
| `const_comp` | Species held constant everywhere in the tube (e.g. `['SO2','O2','H2O']`). Each needs a concentration `modelparams.<X>conc` (from a CSV column, a `<X>ratio` + `<X>flow`, or set in the script) |
| `Init_comp`, `Init_set` | Species fixed at the tube inlet when `Init_set = 'on'` ([4.2](#42-oh-source-point-or-continuous)) |
| `key_spe_for_plot` | Species used for the convergence criterion and progress output (e.g. `'H2SO4'`, `'HOI'`, `'IDHDP'`) |
| `plot_spec` | Species to plot |
| `file_name` | Input CSV name in `Input_files/`; also the name of the output sub-folder |
| `input_file_folder`, `export_file_folder`, `input_mechanism_folder` *(optional)* | Override the default folders |

After the parameter block, `calculate_concentrations(modelparams)` turns the CSV into concentrations. You can add or change concentrations afterwards — Example 4 adds a constant 100 ppb CO background:

```python
calculate_concentrations(modelparams)
CO_conc = 100e-9 * modelparams.p / (1.3806488e-23 * modelparams.TEMP * 1e6)
modelparams.COconc = np.full_like(modelparams.O2conc, CO_conc)
modelparams.const_comp = list(modelparams.const_comp) + ['CO']
modelparams.const_comp_conc = np.transpose([getattr(modelparams, i + 'conc') for i in modelparams.const_comp])
```

All stages with OH > 0 are run (`num_stage = modelparams.OHconc`). To run fewer stages, remove rows from the input CSV.

### 4.5 Model modes

| `model_mode` | Description |
|---|---|
| `'flowtube1'` | Chemistry is integrated together with diffusion/advection at every time step (explicit). Use it for a **point OH source** and small mechanisms. `dt` must resolve the fastest reaction. |
| `'flowtube2'` | Chemistry (stiff solver, parallel over grid cells) and transport are alternated, chemistry over `dt × timesteps` per iteration. Use it for a **continuous OH source** and large mechanisms. |
| `'box'` | 0-D box model: chemistry only, no diffusion or advection. Runs the first CSV row. |
| `'kinetic'` | Transport only, no chemistry (tests the convection–diffusion core). |

### 4.6 Photolysis scaling (SUN)

Mechanisms can contain photolysis reactions written as `SUN * k`, e.g. in `isoprene_reduced_plus_v5_SUNfinal.eqn`:

```
% SUN*5e-5: H2O2 = OH + OH;
```

`SUN` takes the value of `sun_a` from the start script; leaving `sun_a` out (as in Example 3) or setting `sun_a = 0` switches these reactions off. Example 4 uses `sun_a = 2e5`; to scan the photolysis strength, copy the script and change `sun_a` and the input/output name for each case.

### 4.7 Speed vs. accuracy

The examples use coarse grids and loose convergence so that they finish in minutes. Typical settings:

| | Example settings (fast) | Publication settings (example) |
|---|---|---|
| Point source (`flowtube1`, Examples 1–2) | `Rgrid1 = 20`, `Zgrid1 + Zgrid2 = 20`, `dt = 1e-4` (SA) / `1e-3` (HOI), `inter_acc = 1e-3` | `Rgrid1 = 40`, `Zgrid1 = Zgrid2 = 20`, `dt = 5e-5`, `inter_acc = 1e-5` (SA: ~1 min); check with a finer grid |
| Continuous source (`flowtube2`, Example 4) | `Rgrid = 14`, `dt = 2e-3`, `timesteps = 250`, `inter_acc = 1e-3` | `Rgrid = 60`, `Zgrid1 = Zgrid2 = 30`, `Zgridl = 10`, `dt = 5e-5`, `timesteps = 10000`, `inter_acc = 1e-5` (hours to days, use a cluster) |

Grid test for Example 1 (H2SO4, stage 2): R = 20 → 3.46×10⁷, R = 40 → 3.31×10⁷, R = 80 → 3.23×10⁷ cm⁻³ (identical to the original model at R = 80). The fast R = 20 example is therefore ~7 % high and R = 40 ~2.5 % high — use at least R = 40 for calibration factors. Example 4 on the coarse grid agrees with the fine-grid (R = 60) run within ~5 % for OH and HO2 and within ~10–40 % for the oxidation products, so always check grid convergence before using results quantitatively. In `flowtube1` the chemistry is explicit: if the run produces `NaN`, reduce `dt` (fast reactions such as HSO3 + O2 in the SA mechanism need `dt ≈ 1e-4` s).

## 5. Outputs

Results are written to `Export_files/<file_name without .csv>/`. File names encode the run settings, e.g. `R14L16inter_acc2e-05dt2e-03timstep2e+02itx1.9e-10flowtube2_...`:

| File | Content |
|---|---|
| `*_allstage_final_results.csv` | Mean concentration (molecules cm⁻³) at the tube exit of **every species** (columns) for each stage (rows) — the main result |
| `*_Stage_<j>_endtube_profile_output.txt` | Radial profile of all species at the tube end for stage *j* |
| `*_Stage_<j>_1st_tube_profile_output.txt` | Same, at the end of the first tube (two-tube set-ups) |
| `*_Stage_<j>_delta_keyspecforplot.csv` | Convergence history of `key_spe_for_plot`; box mode writes `<file_name>R…box.csv` |
| `*_params.csv` | All parameters used for the run (for reproducibility) |
| `warmstart_stage<j>_*.npz` / `.npy` | Converged 2-D field, used to warm-start a re-run |
| `*checkpoint*.npz` | Periodic checkpoint; a re-run with the same settings **resumes** from it automatically |

**Restart files** (`warmstart_*`, `*checkpoint*`) and **stage carry-over** (stage N+1 starts from the converged field of stage N) are only used for a **continuous OH source**, where fine-grid runs can take hours to days and may need to be resumed. Point-source calibrations are short and always start from a clean tube, so their results never depend on a previous run. Override with `use_restart = True/False` in the start script.

Restart files are named by stage and grid only, **not by the input values** — delete the `*.npz` / `*.npy` files in the output folder after changing the input CSV or concentrations, otherwise a continuous-source run starts from the old field.

## 6. Running on an HPC cluster

Fine grids and large mechanisms are expensive. `run_flowtube_slurm.sh` is a SLURM template; edit the account and environment lines, then:

```bash
cd PANDA520_flowtube
mkdir -p logs
sbatch run_flowtube_slurm.sh Start_SetParam_Isoprene_continuousOH.py
```

The `flowtube2` chemistry step is parallelised over grid cells with Numba — set `NUMBA_NUM_THREADS` to the number of CPUs (done in the template) and keep `OMP/MKL/OPENBLAS_NUM_THREADS=1`. For parameter scans, make one start script per case (different `file_name`) and submit them as a job array. If several jobs share the same folder, give each its own copy of the code (the auto-generated `Funcs/*.py` files are rewritten at the start of each run).

## 7. Tips and troubleshooting

- **`NaN` or exploding concentrations** — reduce `dt`, or coarsen the grid.
- **Point source gives (almost) no product** — check `model_mode = 'flowtube1'`, `Init_comp = ['OH', 'HO2']`, `Init_set = 'on'` and that `Rl` is not set ([4.2](#42-oh-source-point-or-continuous)).
- **No convergence / too slow** — start with a coarse grid and a looser `inter_acc`, then refine. For continuous-source runs, checkpoints and warm starts let you continue interrupted runs: just submit the same script again.
- **`ModuleNotFoundError: Funcs`** — run the scripts from inside `PANDA520_flowtube/` (or add it to `PYTHONPATH`).
- **Open Babel import error** — install it with `conda install -c conda-forge openbabel` (pip wheels are often unavailable).
- **Species not found** — names in `const_comp`, `Init_comp`, `plot_spec`, `key_spe_for_plot` and `Diff_setname` must match the mechanism exactly; species in `const_comp`/`Init_comp` need a concentration.
- **Multiprocessing on macOS/Windows** — keep the `if __name__ == "__main__":` block and `multiprocessing.set_start_method("spawn")` from the examples.

## 8. Acknowledgements

We thank the ACCC Flagship funded by the Academy of Finland grant number 337549, Academy professorship funded by the Academy of Finland (grant no. 302958), Academy of Finland projects no. 346370, 325656, 316114, 314798, 325647, 341349 and 349659. European Research Council (ERC) project ATM-GTP Contract No. 742206. The Arena for the gap analysis of the existing Arctic Science Co-Operations (AASCO) funded by Prince Albert Foundation Contract No 2859. M.K. thanks the Jane and Aatos Erkko Foundation for providing funding. M.K. and X.-C.H thank the Jenny and Antti Wihuri Foundation for providing funding for this research.

## 9. License

MARFORCE-Flowtube is free software released under the [GNU General Public License v3.0](LICENSE). It includes code adapted from [PyCHAM](https://github.com/simonom/PyCHAM) (© 2018–2024 Simon O'Meara, GPL-3.0); those files keep their original copyright headers. You may redistribute and modify it under the terms of the GPL v3; it is distributed WITHOUT ANY WARRANTY.
