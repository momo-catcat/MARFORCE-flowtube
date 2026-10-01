# MARFORCE-Flowtube

MARFORCE-Flowtube simulates **laminar flow tubes with chemistry**: gas flows through a tube with a parabolic velocity profile, species diffuse radially and axially, are lost to the wall and react with each other. The model is mainly used to

- **calibrate chemical-ionisation mass spectrometers (CIMS)** for H2SO4 (SA), HOI and other species: OH is produced by H2O photolysis with a UV lamp, reacts with an excess reactant (SO2, I2, …) and the model predicts the concentration of the product that reaches the instrument;
- simulate **flow-reactor oxidation experiments** (e.g. isoprene chemistry with the Wennberg mechanism), with a continuously illuminated OH-production section and photolysis reactions;
- run the same chemistry as a **0-D box model** for comparison.

The model was published with our manuscript in [Atmospheric Measurement Techniques](https://amt.copernicus.org/articles/16/4461/2023/) — please cite it if you use MARFORCE-Flowtube. Parts of the chemistry parser are adapted from [PyCHAM](https://github.com/simonom/PyCHAM) (Simon O'Meara, GPL-3.0). A Matlab-based model for sulphuric acid calibration only is available at https://github.com/ceciliarighi/ACTRIS_CiGas_condensable_vapors (see [section 7](#7-outlet-concentration-mean-or-flow-weighted) for how its output differs).

Feedback is very welcome — please open an issue.

---

## Contents

1. [Repository layout](#1-repository-layout)
2. [Installation](#2-installation)
3. [Examples — start here](#3-examples--start-here)
4. [How the model works](#4-how-the-model-works)
5. [Setting up your own simulation](#5-setting-up-your-own-simulation)
6. [Outputs](#6-outputs)
7. [Outlet concentration: mean or flow-weighted](#7-outlet-concentration-mean-or-flow-weighted)
8. [Long runs, restarts and HPC](#8-long-runs-restarts-and-hpc)
9. [Troubleshooting](#9-troubleshooting)
10. [Citation, acknowledgements and license](#10-citation-acknowledgements-and-license)

---

## 1. Repository layout

```
MARFORCE-flowtube/
├── README.md                    # this manual
├── LICENSE                      # GNU GPL v3
├── environment.yml              # conda environment (recommended)
├── requirements.txt             # pip packages (Open Babel via conda)
└── PANDA520_flowtube/           # run everything from this folder
    ├── Start_SetParam_SA_calibration_direct.py   # Example 1a: SA calibration, point OH, direct chemistry
    ├── Start_SetParam_SA_calibration_ODE.py      # Example 1b: same, ODE-solver chemistry
    ├── Start_SetParam_HOI_calibration.py         # Example 2:  HOI calibration, point OH, two tube sections
    ├── Start_SetParam_Isoprene_box.py            # Example 3:  isoprene chemistry, box model
    ├── Start_SetParam_Isoprene_continuousOH.py   # Example 4:  isoprene, continuous OH + SUN photolysis
    ├── Start_SetParam_kinetic_test.py            # Example 5:  transport only (no chemistry) vs Gormley-Kennedy
    ├── run_flowtube_slurm.sh                     # SLURM job template
    ├── Input_files/          # input CSVs: one row = one experiment stage
    ├── input_mechanism/      # chemical mechanisms (SA, HOI, Isoprene/Wennberg, Isoprene/MCM)
    ├── photofiles/           # MCM v3.2 absorption cross-sections / quantum yields, lamp spectrum
    ├── Funcs/                # model source code (see section 4.3)
    ├── kinetics/             # diffusion coefficients and helper kinetics
    └── Export_files/         # results (one sub-folder per input CSV)
```

## 2. Installation

1. **Get the code**

   ```bash
   git clone https://github.com/momo-catcat/MARFORCE-flowtube.git
   cd MARFORCE-flowtube
   ```

2. **Create the Python environment** (conda is recommended because Open Babel is easiest to install from conda-forge):

   ```bash
   conda env create -f environment.yml
   conda activate marforce
   ```

   or, with an existing Python 3.9–3.10:

   ```bash
   conda install -c conda-forge openbabel=3.1.1
   pip install -r requirements.txt
   ```

   | Package | Version tested | Used for |
   |---|---|---|
   | numpy | 2.0.2 | arrays |
   | scipy | 1.13.1 | sparse Jacobians, interpolation, constants |
   | pandas | 2.2.3 | input CSVs, result tables |
   | numba | 0.60.0 | compiled, parallel transport and chemistry kernels |
   | numbalsoda | 0.3.4 | LSODA stiff ODE solver called from Numba (`flowtube2`, `box`) |
   | matplotlib | 3.9.4 | concentration plots |
   | molmass | 2024.10.25 | molar masses / diffusion coefficients from formulas |
   | xmltodict | 0.14.2 | reading `chemical_species_custom.xml` |
   | openbabel | 3.1.1 | SMILES → molecular properties |

3. **Check the installation**

   ```bash
   python -c "import numpy, scipy, pandas, numba, numbalsoda, matplotlib, molmass, xmltodict; from openbabel import pybel; print('OK')"
   ```

## 3. Examples — start here

Run every script from inside `PANDA520_flowtube/`:

```bash
cd PANDA520_flowtube
python Start_SetParam_SA_calibration_direct.py
```

Each flow-tube example simulates **two experiment stages** (two rows of its input CSV). Plots are updated while the model runs; set `export MPLBACKEND=Agg` on a machine without a display. Results are written to `Export_files/<input CSV name>/` ([section 6](#6-outputs)).

| # | Script | Set-up | OH source | Solver | Grid (R × Z) | Run time* |
|---|---|---|---|---|---|---|
| 1a | `Start_SetParam_SA_calibration_direct.py` | SA calibration, one 3/4" tube (R = 0.78 cm, 26 cm) | point | `flowtube1` (direct) | 80 × 40 | ~1 min |
| 1b | `Start_SetParam_SA_calibration_ODE.py` | same as 1a | point | `flowtube2` (ODE) | 80 × 40 | ~0.7 min |
| 2 | `Start_SetParam_HOI_calibration.py` | HOI calibration, 0.78 cm (41 cm) + 1.04 cm (58.5 cm) sections | point | `flowtube1` | 80 × (20 + 20) | ~3 min |
| 3 | `Start_SetParam_Isoprene_box.py` | isoprene oxidation, chamber-like box | held constant | `box` | – | ~0.5 min |
| 4 | `Start_SetParam_Isoprene_continuousOH.py` | isoprene oxidation, lamp section + two tubes with Y-piece, SUN photolysis, 100 ppb CO | continuous | `flowtube2` | 14 × (6 + 10 + 10) | ~4 min |
| 5 | `Start_SetParam_kinetic_test.py` | transport test: H2SO4 through a 1 m tube without chemistry, compared with Gormley–Kennedy | – | `kinetic` | 80 × 40 | ~0.5 min |

\*Laptop with 10 CPU cores.

**Results of the examples** (outlet concentrations, molecules cm⁻³, `final_output_method = 'mean'`):

| Example | Stage 1 | Stage 2 |
|---|---|---|
| 1a SA, direct | H2SO4 = 6.74×10⁷, OH = 2.58×10⁶ | H2SO4 = 3.23×10⁷, OH = 1.26×10⁶ |
| 1b SA, ODE | H2SO4 = 6.73×10⁷ | H2SO4 = 3.23×10⁷ |
| 2 HOI | HOI = 6.76×10⁶ | HOI = 8.57×10⁶ |
| 3 Isoprene box | IDHDP = 2.80×10⁹ | – (box mode runs one stage) |
| 4 Isoprene continuous OH | IDHDP = 2.33×10⁶, OH = 5.19×10⁸ | IDHDP = 1.28×10⁶, OH = 4.07×10⁸ |
| 5 Transport test (`'weighted'`) | H2SO4 penetration 0.657 (Gormley–Kennedy 0.645, +1.9 %) at 22.5 slpm | 0.458 (0.449, +2.0 %) at 10 slpm |

### 3.1 Which SA calibration example should I use? (1a vs 1b)

Both examples simulate the identical calibration; they differ only in how the chemistry is integrated ([section 4.2](#42-numerical-method)):

| | **1a — direct (`flowtube1`)** | **1b — ODE solver (`flowtube2`)** |
|---|---|---|
| Chemistry | computed explicitly inside every transport time step (the original MARFORCE method) | stiff ODE solver (LSODA), parallel over grid cells, alternating with transport every 1 ms |
| Time step | limited by the fastest reaction (HSO3 + O2 → `dt = 1e-4` s) | limited by transport stability only (`dt = 2.5e-4` s) |
| H2SO4 stage 1 / stage 2 | 6.74×10⁷ / 3.23×10⁷ | 6.73×10⁷ / 3.23×10⁷ |
| Run time (R80 × Z40, 2 stages) | 1.1 min | 0.7 min |
| Best for | small mechanisms; simplest and closest to the published model | large or stiff mechanisms (many species, fast reactions), many grid cells |

The two methods agree within **0.2 %**, and the result reproduces the original MARFORCE code exactly (H2SO4 = 6.735×10⁷ cm⁻³ for the same input). The ODE result does not depend on the 1-ms alternation interval (0.5, 1 and 2 ms give 6.74, 6.73 and 6.73×10⁷). **Use 1a** for standard SA/HOI calibrations; **use 1b** when the mechanism is too stiff for small explicit time steps.

**Grid resolution** (Example 1a, H2SO4 stage 2): R20 → 3.46×10⁷, R40 → 3.31×10⁷, R80 → 3.23×10⁷. A coarse grid over-estimates H2SO4 by up to ~7 %; use R ≥ 80 (the original resolution) for calibration factors.

### 3.2 Checking the transport core (Example 5)

`model_mode = 'kinetic'` switches the chemistry off: species are only carried by the laminar flow, diffuse and are lost to the wall. Example 5 puts 10⁸ cm⁻³ H2SO4 at the inlet of a 1 m tube and compares the flow-weighted fraction that leaves the tube with the analytical Gormley & Kennedy (1949) penetration `P(μ)`, `μ = πDL/Q` (the formula of the original Matlab calibrator). The model agrees within 2 % at both flows; refining the axial grid changes this by < 0.3 %, the remaining difference comes from the approximations of the analytical formula (no axial diffusion). Use this mode to check grid and `dt` of a new set-up before adding chemistry: copy Example 5, change geometry, flows and `Diff_set`, and compare.

## 4. How the model works

### 4.1 Physics

For every species the model solves the axisymmetric convection–diffusion–reaction equation

```
∂c/∂t = D (∂²c/∂r² + (1/r) ∂c/∂r + ∂²c/∂z²)  −  u(r) ∂c/∂z  +  P(c) − L(c)
```

- **Laminar flow:** `u(r) = 2Q/(πR²) · (1 − r²/R²)` (Poiseuille profile; `Q` volume flow, `R` tube radius).
- **Diffusion:** `D` from `Diff_set` in the start script, otherwise estimated from the molecular formula with Fuller's method (`kinetics/diff_coef.py`).
- **Chemistry:** `P − L` from the chemical mechanism (rate coefficients may depend on T, p, M, O2, H2O, photolysis scaling `SUN`, lamp term `itProd`).
- **Boundary conditions:** at the **wall** the concentration of every reactive species is 0 (the wall is a perfect sink — diffusion-limited wall loss); **constant species** (`const_comp`) are held fixed everywhere; at the **inlet** the `Init_comp` species are fixed (when `Init_set = 'on'`), all other reactive species enter with zero concentration; the **outlet** is open.
- **OH source:**
  - *point* — OH and HO2 are created in a short lamp region at the inlet: `[OH]₀ = [HO2]₀ = It · σ(H2O) · [H2O]`, with `It = Itx · Qx / Q` (lamp It product scaled to the actual flow) and σ(H2O, 185 nm) = 7.22×10⁻²⁰ cm²;
  - *continuous* — OH is produced along an illuminated section of length `Ll` by the mechanism reaction `H2O = OH + H` with first-order rate `Itx` (s⁻¹).
- **Several tubes:** a 2nd tube (radius `R2`, length `L2`) is simulated after the 1st when its radius or its flow differs (Y-piece). The outlet profile of tube 1 is mapped ring by ring onto the inlet of tube 2 and diluted by `Q1/Q2`.

The model integrates in time until the concentrations no longer change, i.e. it computes the **steady state** of the tube. The **outlet concentration** is then averaged over the cross-section ([section 7](#7-outlet-concentration-mean-or-flow-weighted)).

### 4.2 Numerical method

- **Grid:** `Rgrid` points across the full diameter × `Zgrid` points along the tube (finite differences). The solution is symmetric about the axis.
- **Iterations:** the model repeatedly advances the solution by `timesteps` steps of length `dt` (one *iteration* = `dt × timesteps` seconds). After each iteration it compares `key_spe_for_plot` at the outlet with the previous iteration; the run has converged when the relative change is below `inter_acc` (and at least `fix_timstep` iterations were done).
- **Stages:** every row of the input CSV is a separate steady-state calculation (e.g. different H2O flows).
- **Solvers (`model_mode`):**

  | `model_mode` | How chemistry and transport are combined | Time-step limit |
  |---|---|---|
  | `'flowtube1'` | **direct**: diffusion, advection and chemistry are advanced together with explicit Euler steps | `dt` < lifetime of the fastest-reacting species and the transport stability limit |
  | `'flowtube2'` | **operator splitting**: every iteration first integrates the chemistry in every grid cell for `dt × timesteps` with a stiff ODE solver (LSODA, Numba-parallel), then transports for `timesteps` steps of `dt` | transport stability only; `dt × timesteps` must be short compared with the residence time |
  | `'box'` | 0-D: chemistry (stiff solver), dilution and wall loss only, until steady state | none |
  | `'kinetic'` | transport only: the chemistry step is skipped, species set at the inlet are only advected, diffused and lost to the wall (Example 5) | transport stability |

- **Transport stability (explicit scheme):** roughly `dt < dr²/(4D)` (dr = 2R/(Rgrid−1)) and `dt < dx/u_max`. Finer grids need smaller `dt`; an unstable `dt` produces `NaN` and the run stops with a message.

### 4.3 Program flow and code map

```
Start_SetParam_*.py            all parameters (SimulationParams)
 └─ calculate_concentrations   Funcs/Calcu_by_flow.py   input CSV → concentrations of every stage,
 │                                                       [OH]0, tube configuration, file paths
 └─ Run_flowtube               Funcs/Run_flowtube.py    loop over stages, write result files
     └─ cmd_calib5             Funcs/cmd_calib5.py      one stage:
         ├─ eqn_pars.extr_mech   parse the mechanism, write Funcs/rate_coeffs.py, ode_solv.py, ... (auto-generated)
         ├─ grid_para            grid (Funcs/grid_parameters.py), diffusion coefficients (get_diff_and_u.py)
         ├─ model_onetube / model_twotubes / model_box      iterate to steady state
         │     ├─ ode_solv_batch / ode_solv_numba_batch      chemistry with LSODA (flowtube2, box)
         │     ├─ odesolve3.odesolve                         transport (+ chemistry in flowtube1)
         │     └─ set_boundlay_for_ode, cal_const_comp_conc  boundary conditions, constant species
         └─ meanconc_cal          outlet concentration ('mean' or 'weighted')
```

| File | Purpose |
|---|---|
| `Funcs/Calcu_by_flow.py` | Flows, bottle mixing ratios and H2O saturation → concentrations; point-source [OH]₀; tube type |
| `Funcs/Run_flowtube.py` | Runs all stages and writes the result tables |
| `Funcs/cmd_calib5.py` | Sets up and runs one stage; chooses one-tube / two-tube / box model; restart files |
| `Funcs/eqn_pars.py`, `eqn_interr.py`, `sch_interr.py`, `xml_interr.py` | Mechanism parser (adapted from PyCHAM) |
| `Funcs/write_*.py` | Write the auto-generated `rate_coeffs.py`, `ode_solv.py`, `dydt_rec.py`, `hyst_eq.py` (do not edit those by hand) |
| `Funcs/model_onetube.py`, `model_twotubes.py`, `model_box.py` | Iteration loops, convergence, checkpoints |
| `Funcs/odesolve3.py` | Finite-difference transport (and explicit chemistry for `flowtube1`) |
| `Funcs/ode_solv_batch.py`, `ode_solv_numba_batch.py`, `ode_worker.py` | Parallel stiff chemistry solver |
| `Funcs/meanconc_cal.py` | Outlet concentration (mean / flow-weighted) |
| `Funcs/Wennberg_rec_funcs.py`, `itprod.py`, `kclust.py` | Rate-coefficient helper functions usable in mechanisms |
| `Funcs/read_output_profile_file.py` | Read the saved 2-D concentration fields ([section 6](#6-outputs)) |
| `kinetics/diff_coef.py` | Diffusion coefficients (Fuller's method) |

## 5. Setting up your own simulation

1. **Copy the closest example** and give it a new name, e.g. `cp Start_SetParam_SA_calibration_direct.py Start_SetParam_mycal.py`.
2. **Write the input CSV** for your experiment in `Input_files/` ([5.1](#51-input-csv-and-tube-configurations)) and set `file_name` to it.
3. **Set the geometry and the OH source** ([5.2](#52-oh-source-point-or-continuous)): radii, lengths, `Itx`, `Qx`.
4. **Choose the mechanism** ([5.3](#53-chemical-mechanism)) and the species treatment: `const_comp`, `Init_comp`, `key_spe_for_plot`.
5. **Choose grid, `dt` and convergence** ([5.6](#56-choosing-grid-time-step-and-convergence)). Start coarse, then refine until the result stops changing.
6. **Run** `python Start_SetParam_mycal.py` from `PANDA520_flowtube/` and read the results in `Export_files/<CSV name>/`.

### 5.1 Input CSV and tube configurations

One row = one **experiment stage**; only include the stages you want to run. Columns are recognised by their names:

| Column | Meaning |
|---|---|
| `<X>flow` | Flow of gas X in **sccm**. A trailing `1`/`2` means tube 1 / flow added at the Y-piece (e.g. `H2Oflow1`, `N2flow2`) |
| `<X>conc` | Measured or prescribed concentration of X (molecules cm⁻³), e.g. `I2conc`, `ISOP1OH2OOHconc`, `H2Oconc`. Available as `modelparams.<X>conc` for `const_comp` / `Init_comp` |
| `Q` or `Q1`, `Q2` | Total flow (sccm) in the tube / in tube 1 and tube 2. `Q2` is the total flow in tube 2, not the Y-piece flow |
| `T` *(optional)* | Temperature per stage (K) |
| `time` *(optional)* | Run time per stage (s); raises the minimum number of iterations |
| others | ignored (e.g. `UVC`, timestamps) |

How concentrations are calculated from flows (`Calcu_by_flow.py`):

- gas from a bottle: `[X] = Xflow × Xratio / Q × p/(k_B T)`, with `<X>ratio` the mole fraction in the bottle given in the start script (e.g. `SO2ratio = 5e-3`);
- O2: `O2flow × O2ratio / Q × p/(k_B T)` (O2flow is synthetic air);
- H2O: `H2Oflow / Q × p_sat(T)/(k_B T)` (saturated humidifier). If an `H2Oconc` column is given it is used instead **when the H2O flow is ≥ 1000 sccm**, also for [OH]₀;
- `outflowLocation = 'before'` uses `Q` for the dilution, `'after'` uses the sum of the injected flows.

**Tube configurations** (all work with point and continuous OH sources):

| Configuration | What to give | Model used |
|---|---|---|
| One tube | single-flow CSV (`H2Oflow`, `N2flow`, `Q`); `R1`, `L1`; no `L2` | one tube |
| Two sections with **different diameter**, same flow — e.g. 3/4" (R = 0.78 cm) followed by 1" (R = 1.2 cm) | single-flow CSV; `R1`, `L1`, `R2`, `L2` | two tubes (Example 2) |
| Two tubes with a **Y-piece** adding flow | two-flow CSV (`H2Oflow1`, `N2flow1`, `Q1`, `Q2`, `…flow2`); `R1`, `L1`, `R2`, `L2` | two tubes with dilution `Q1/Q2` (Example 4) |

Keep `Rgrid2 = Rgrid1`. For a one-tube set-up the axial grid is `Zgrid1 + Zgrid2` points over `L1`; for two tubes tube 1 has `Zgrid1` and tube 2 `Zgrid2` points. Note that in fully developed laminar flow the diffusional wall loss depends on `D·L/Q` and not on the radius, so changing only the diameter mainly matters through the longer residence time in a wider tube (more reaction time).

### 5.2 OH source: point or continuous

The model decides from `Rl`: **`Rl > 0` → continuous source, otherwise point source.**

| | **Point source** (Examples 1, 2) | **Continuous source** (Example 4) |
|---|---|---|
| Physical picture | OH/HO2 formed in a short lamp region at the inlet, then react downstream (classic calibration) | OH produced along an illuminated section while the chemistry proceeds |
| Geometry | `R1`, `L1` (+ `R2`, `L2`); **no `Rl`, `Ll`** | `Rl`, `Ll`, `Rgridl`, `Zgridl` + `R1`, `L1`, `R2`, `L2` |
| `Itx` | lamp *It* product at `Qx` (photons cm⁻², e.g. `4.84e10`) | first-order rate coefficient of `H2O = OH + H` (s⁻¹, e.g. `1.9e-10`) |
| Mechanism | – | must contain `% itProd(modelparams) : H2O = OH + H ;` |
| `Init_comp`, `Init_set` | `['OH', 'HO2']`, `'on'` | OH not in `Init_comp`; `'on'` only for other inlet species (e.g. `['ISOP1OH2OOH']`) |
| `model_mode` | `'flowtube1'` (or `'flowtube2'` with `dt × timesteps` ≈ 1 ms, Example 1b) | `'flowtube2'` |

### 5.3 Chemical mechanism

Put the mechanism in its own folder under `input_mechanism/` and set `folder_mechaism` and `sch_name`.

| Format | `flag_mech` | `tsv_file` |
|---|---|---|
| MCM export (`.fac`, from the [MCM website](https://mcm.york.ac.uk/MCM)) | `'1'` | the MCM species export `mcm_export_species.tsv`; `chemical_species_custom.xml` is generated automatically |
| Other mechanisms (`.eqn`, `.txt`) | `'0'` | a species table (e.g. `SMILE_formula.csv`); provide `chemical_species_custom.xml` (species → SMILES) in the same folder |

File format (see `input_mechanism/HOI/HOI_cali_chem.txt` for a short, complete example):

```
* lines starting with * are comments ;
KMT06 = 1 + (1.40D-21*EXP(2200/TEMP)*H2O) ;          <- generic rate coefficient (name < 10 characters)
% 2.1D-10 : I2 + OH = HOI + I ; # reference           <- reaction:  % rate coefficient : reactants = products ;
```

Rate expressions may use Fortran or Python notation, `EXP/exp, LOG, LOG10, dsqrt, dabs, numpy.*`, the variables `TEMP` (K), `p` (Pa), `M`, `N2`, `O2`, `H2O` (molecules cm⁻³), the helper functions `TROE, TUN, ALK, NIT, EPO, ISO1, ISO2, KCO` (Wennberg), `kDimer/kTrimer` (H2SO4 clustering), `itProd(modelparams)` (continuous OH source), `SUN` ([5.5](#55-photolysis-scaling-sun)) and MCM photolysis rates `J(n)`. Do not leave blank lines in the mechanism file.

### 5.4 Parameter reference

All parameters are set in the `SimulationParams(...)` call of the start script (numbers are converted to `float32`, grid sizes to `int`). Every example script documents its parameters line by line.

| Group | Parameter | Description |
|---|---|---|
| Conditions | `p`, `TEMP` | Pressure (Pa), temperature (K) |
| Geometry | `R1`, `L1`, `R2`, `L2` | Inner **radius** (cm) and length (cm) of tube 1 / tube 2 |
| | `Rl`, `Ll` | Radius and length of the illuminated section (continuous source only) |
| OH source | `Itx`, `Qx` | Lamp parameter and the flow (slpm) at which it was determined ([5.2](#52-oh-source-point-or-continuous)) |
| | `sun_a` | Value of `SUN` in the mechanism; default 0 |
| Flows | `sampleflow` | CIMS inlet flow (slpm) |
| | `outflowLocation` | `'before'`/`'after'` the gas injection |
| | `O2ratio`, `<X>ratio` | O2 fraction of synthetic air; mole fraction of X in its gas bottle |
| Grid | `Rgrid1`, `Zgrid1`, `Rgrid2`, `Zgrid2` | Radial (across the diameter) / axial points of tube 1 / tube 2 |
| | `Rgridl`, `Zgridl` | Grid of the illuminated section |
| Numerics | `model_mode` | `'flowtube1'`, `'flowtube2'`, `'box'`, `'kinetic'` ([4.2](#42-numerical-method)) |
| | `dt`, `timesteps` | Time step (s) and steps per iteration |
| | `inter_acc`, `fix_timstep` | Convergence threshold and minimum number of iterations. In box mode the threshold is `inter_acc × timesteps / 10⁴` |
| | `num_plot` | Print/plot every *n* iterations |
| | `Diff_setname`, `Diff_set` | User diffusion coefficients (cm² s⁻¹) |
| | `use_restart` | Restart files ([section 8](#8-long-runs-restarts-and-hpc)); default on for continuous, off for point sources |
| Chemistry | `folder_mechaism`, `sch_name`, `tsv_file`, `flag_mech` | Mechanism ([5.3](#53-chemical-mechanism)) |
| | `const_comp` | Species held constant everywhere; each needs a concentration (`<X>conc` column, `<X>ratio` + `<X>flow`, or set in the script) |
| | `Init_comp`, `Init_set` | Species fixed at the inlet when `Init_set = 'on'` |
| | `key_spe_for_plot` | Species for the convergence test |
| | `plot_spec` | Species to plot |
| Output | `file_name` | Input CSV in `Input_files/`; also the output folder name |
| | `final_output_method` | `'mean'` (default) or `'weighted'` ([section 7](#7-outlet-concentration-mean-or-flow-weighted)) |
| | `input_file_folder`, `export_file_folder`, `input_mechanism_folder` | Optional folder overrides |
| Box model | `dil_fac_now` | Dilution rate (s⁻¹) (or `chamber_volume` [L] + `total_flow` [slpm]) |
| | `wall_loss_set`, `wall_loss_custom` | Wall loss: 0 = off, 1 = estimated, x = x s⁻¹; per-species values (s⁻¹) |
| Internal | `tot_time`, `save_step`, `pars_skip` | Keep the example values (`pars_skip = 0` parses the mechanism) |

Concentrations can be added or changed after `calculate_concentrations(modelparams)`; Example 4 adds a constant CO background this way.

### 5.5 Photolysis scaling (SUN)

Reactions written as `SUN*k` in the mechanism, e.g. `% SUN*5e-5: H2O2 = OH + OH;`, use `SUN = sun_a` from the start script. Leave `sun_a` out or set it to 0 to switch them off (Example 3); Example 4 uses `sun_a = 2e5`.

### 5.6 Choosing grid, time step and convergence

1. **Grid:** start coarse (R 20–40) to set up the case, then refine until the result changes by less than your uncertainty. For calibration factors use R ≥ 80 (Example 1a: R20 is 7 % and R40 2.5 % above R80).
2. **`dt`:**
   - `flowtube1`: below the lifetime of the fastest-reacting species (SA: HSO3, ~2×10⁻⁴ s → `dt = 1e-4`) **and** below the transport limit of the grid (HOI on R80: `dt = 1e-4`).
   - `flowtube2`: transport limit only. Keep `dt × timesteps` short compared with the residence time for a point source (Example 1b: 1 ms).
   - Check: halving `dt` must not change the result. `NaN` means `dt` is too large.
3. **Convergence:** `inter_acc` is the relative change of `key_spe_for_plot` between iterations. Use 1e-4 or smaller for calibrations; one iteration (`dt × timesteps`) should be at least about one residence time for `flowtube1`.

Slow, high-resolution runs of the continuous-source case (publication settings: `Rgrid = 60`, `Zgrid1 = Zgrid2 = 30`, `Zgridl = 10`, `dt = 5e-5`, `timesteps = 10000`, `inter_acc = 1e-5`) take hours to days — use a cluster ([section 8](#8-long-runs-restarts-and-hpc)). On the coarse grid of Example 4 OH and HO2 agree with that fine grid within ~5 %, oxidation products within ~10–40 %.

## 6. Outputs

Results go to `Export_files/<file_name without .csv>/`; file names encode the settings (grid, convergence, `dt`, `timesteps`, `Itx`, mode):

| File | Content |
|---|---|
| `*_allstage_final_results.csv` | **Main result:** outlet concentration (molecules cm⁻³) of every species (columns) for every stage (rows) |
| `*_Stage_<j>_delta_keyspecforplot.csv` | Convergence history: time, relative change and outlet concentration of all species per iteration |
| `*_Stage_<j>_endtube_profile_output.txt` | Full 2-D concentration field (all grid cells, all species) at convergence (last tube) |
| `*_Stage_<j>_1st_tube_profile_output.txt` | Same for tube 1 (two-tube set-ups) |
| `*_params.csv` | All parameters of the run |
| `warmstart_*`, `*checkpoint*` | Restart files (continuous source only, [section 8](#8-long-runs-restarts-and-hpc)) |

Reading the results in Python (run from `PANDA520_flowtube/`):

```python
import glob, pandas as pd
from types import SimpleNamespace
from Funcs.read_output_profile_file import read_output_profile

folder = 'Export_files/SA_calibration/'
res = pd.read_csv(glob.glob(folder + '*allstage_final_results.csv')[0], index_col=0)
print(res[['OH', 'H2SO4']])                       # outlet concentrations per stage

# 2-D field of stage 0: array [radial index, axial index, species]
prm = pd.read_csv(glob.glob(folder + '*params.csv')[0]).set_index('Parameter')['Value']
grid = SimpleNamespace(Rgrid=int(prm['Rgrid']), Zgrid=int(prm['Zgrid']), comp_num=len(res.columns))
c = read_output_profile(glob.glob(folder + '*Stage_0_endtube_profile_output.txt')[0], grid)
h2so4_outlet_profile = c[:, -1, list(res.columns).index('H2SO4')]
```

## 7. Outlet concentration: mean or flow-weighted

The outlet profile is not uniform (low at the wall, high in the centre), so it has to be reduced to one number. `final_output_method` selects how:

- **`'mean'` (default)** — area-weighted average over the outlet cross-section. Each ring's concentration is weighted by its area.
- **`'weighted'`** — flow-weighted average: each ring is weighted by its area **and** its local velocity, i.e. the amount of the species leaving the tube per unit time divided by the flow. This is what the Matlab-based calibration model computes.

Our model defaults to the average concentration because the flow at the end of the tube is not laminar: only ~0.8 slpm is sampled into the mass spectrometer and the rest (typically 10–20 slpm, depending on the inlet system) goes to the exhaust, so the gas is mixed before it is sampled. The models assume perfectly laminar flow to the end, which may introduce errors in the quantification; the area average constrains the mass balance of the species before it enters the instrument. **Users can choose either method — both are supported.** For Example 1a (stage 1) the two methods give H2SO4 = 6.74×10⁷ (`'mean'`) and 8.43×10⁷ cm⁻³ (`'weighted'`, +25 %) — both identical to the original MARFORCE code; the methods typically differ by 10–25 %.

## 8. Long runs, restarts and HPC

**Restart files** are written only for a **continuous OH source**, whose fine-grid runs can take hours to days:

- *checkpoints* (`*checkpoint*.npz`) every 30 min and when the job is stopped by SLURM — start the same script again and the run continues;
- *warm start* (`warmstart_stage<j>_*`) — the converged field of a stage, used as the starting point when the run is repeated;
- *stage carry-over* — stage N+1 starts from the converged field of stage N.

Point-source calibrations are short and always start from a clean tube, so their results never depend on earlier runs. Override with `use_restart = True/False`. Restart files are named by stage and grid only, **not by the input values** — delete them after changing the input CSV or concentrations.

**HPC:** `run_flowtube_slurm.sh` is a SLURM template (edit account and environment lines):

```bash
cd PANDA520_flowtube
mkdir -p logs
sbatch run_flowtube_slurm.sh Start_SetParam_Isoprene_continuousOH.py
```

Set `NUMBA_NUM_THREADS` to the number of CPUs (done in the template) and keep `OMP/MKL/OPENBLAS_NUM_THREADS=1`. For parameter scans make one start script per case (different `file_name`) and submit a job array. Give every simultaneous job its own copy of the code, because the auto-generated `Funcs/*.py` files are rewritten at the start of each run.

## 9. Troubleshooting

| Problem | Solution |
|---|---|
| `NaN detected — aborting this run` | `dt` too large for the grid or the chemistry: halve `dt` ([5.6](#56-choosing-grid-time-step-and-convergence)) |
| Point source gives (almost) no product | `Init_comp = ['OH', 'HO2']`, `Init_set = 'on'`, no `Rl`; with `flowtube2` keep `dt × timesteps` ≈ 1 ms |
| Result changes when re-running a continuous-source case | old restart files: delete `warmstart_*` / `*checkpoint*` in the output folder |
| `ModuleNotFoundError: Funcs` | run from inside `PANDA520_flowtube/` |
| Open Babel import error | `conda install -c conda-forge openbabel=3.1.1` |
| Species not found / `KeyError` | names in `const_comp`, `Init_comp`, `plot_spec`, `key_spe_for_plot`, `Diff_setname` must match the mechanism; `const_comp`/`Init_comp` species need a concentration |
| `IndexError` while reading the mechanism | blank line or old `{1.} A = B : k ;` format in the mechanism file — use `% k : A = B ;` |
| Multiprocessing errors on macOS/Windows | keep the `if __name__ == "__main__":` block with `multiprocessing.set_start_method("spawn")` |

## 10. Citation, acknowledgements and license

Please cite the [AMT article](https://amt.copernicus.org/articles/16/4461/2023/) when using MARFORCE-Flowtube.

We thank the ACCC Flagship funded by the Academy of Finland grant number 337549, Academy professorship funded by the Academy of Finland (grant no. 302958), Academy of Finland projects no. 346370, 325656, 316114, 314798, 325647, 341349 and 349659. European Research Council (ERC) project ATM-GTP Contract No. 742206. The Arena for the gap analysis of the existing Arctic Science Co-Operations (AASCO) funded by Prince Albert Foundation Contract No 2859. M.K. thanks the Jane and Aatos Erkko Foundation for providing funding. M.K. and X.-C.H thank the Jenny and Antti Wihuri Foundation for providing funding for this research.

**AI assistance.** Version 2 of MARFORCE-Flowtube was developed with substantial help from [Claude Code](https://claude.com/claude-code), an AI coding assistant by Anthropic, working under the direction of the authors throughout 2026. Claude Code contributed to:

- **numerical development** — vectorised reaction rates and Jacobians, the parallel stiff chemistry solver (Numba + LSODA), multiprocessing, convergence acceleration, and checkpoint / warm-start / restart handling for long HPC runs;
- **debugging and testing** — e.g. finding an error in the diffusion coefficients of H2SO4 clusters (`kinetics/diff_coef.py`), the stage inlet and two-tube fixes, the transport test against Gormley & Kennedy, and comparisons with the original MARFORCE code;
- **calibration and sensitivity studies** — separate dimer/trimer clustering rates (`kDimer`, `kTrimer`), parameter sweeps, comparison scripts and SLURM job scripts;
- **the final release** — code clean-up, example scripts, requirement files and this README.

The scientific design, chemical mechanisms, physical assumptions, parameter choices and the interpretation of results are the authors'. Results were checked by running the examples, against the original MARFORCE code ([section 3.1](#31-which-sa-calibration-example-should-i-use-1a-vs-1b)), against the analytical Gormley–Kennedy solution ([section 3.2](#32-checking-the-transport-core-example-5)) and against measurements. Commits made with this assistance are marked "Co-Authored-By: Claude" in the git history.

MARFORCE-Flowtube is free software released under the [GNU General Public License v3.0](LICENSE). It includes code adapted from [PyCHAM](https://github.com/simonom/PyCHAM) (© 2018–2024 Simon O'Meara, GPL-3.0); those files keep their original copyright headers. It is distributed WITHOUT ANY WARRANTY.
