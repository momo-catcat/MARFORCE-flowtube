# PANDA520-flowtube
This is a project to develop a flow tube chemistry module that include 2D space and chemistry processes

Shall we include the instruction for SA calibration experiments

## Set parameters (e.g., Start_SetParam_SA_example.py)

`sampleflow` : inlet flow of CIMs, unit lpm \
`R1` : inner diameter for the first tube, unit cm \
`L1` : length for the first tube, unit cm \
`R2` : inner diameter for the second tube if there is any, = 0 for `flag_tube = 1` (this will be introduced in detail at **Calculation of concentrations by flows**), unit cm \
`L2` : length for the second tube, unit cm \
`Itx` : It product at `Qx` flow rate \
`SO2ratio` : SO2 ratio of the gas bottle, unit ppm \
`outflowLocation` : outflow (exhaust) tube located 'before' or 'after' injecting synthetic air, water vapor and SO2 \
`input_file_folder` : the folder where the csv file containing input information (such as flow rates of H2O, SO2 and etc.) locates; default folder is */../Input_files/* if this parameter is missing \
`export_file_folder` : the folder where the csv and txt files containing export infomration (such as H2SO4 concentration) locate; default folder is */../Export_files* if this paramter is missing \
`file_name` : the file you store all the data at different stage including flow rates of N2, O2, SO2, H2O and total flow Q and temperature T if recorded in the folder input files \
`p` : pressure, unit Pa \
`T` : temperature, unit K \
`fullOrSimpleModel` : 'simple' means Gormley & Kennedy approximation, while 'full' means flow model (much slower) \
`O2ratio` : O2 ratio in synthetic air \
`Zgrid` : number of grids in direction of tube length \
`Rgrid` : number of grids in direction of radius \
`dt` : differential time interval (usually 1e-4, but try small number if the values are too large out of range and showing 'Nan' during calculation) \
`model_mode` : use 'normal' if you don't know what this is for. 'kinetic' mode refers to running the model without chemistry module to test the kinetic core. \
`Diff_setname` :  diffusion for the species that you want to define by yourself, otherwise it will be calculated automatically based on the elements it contains \
`Diff_set` : add the value according to the `Diff_setname` \
`sch_name` : chemical scheme file name store in the input_mechanism folder ("SO2_SA.txt" for SA calibration and "HOI_cali_chem.txt" for HOI calibration) \
`const_comp` : species you think they should have constant concentrations in the whole tube (for SA calibration, `const_comp=['SO2','O2','H2O']`) \
`Init_comp` : species you think they should have initial concentration in the first grid of tube (for SA calibration, `Init_comp=['OH','HO2']`) \
`key_spe_for_plot` : key species as criterion to stop the loop (for SA calibration, `key_spe_for_plot='H2SO4'`) \
`plot_spec` : species that you want to plot
`flag_tube` : usually this is determined in the concentration calculation function (**def *calculate_concs***), but `flag_tube=4` needs to be set manually here if needed. \
`Init_comp_conc` : set the initial concentrations for species that you already set in the paras (in the case of SA calibration, `OHconc` from next section is used for both 'OH' and 'HO2') \
`const_comp_conc` : set the constant concentrations for species you already set, but also the specific calculation is from next section. \
`num_stage` : default is `num_stage= paras['OHconc']`, meaning the number of stages refer to OHconc since we want to check the stage without light on, but this can also be set by user as a integer number N (in this case, it means the first N stages are calculated)


## Calculation of concentrations by flows (e.g., Calcu_by_flow_SA.py)

### **Function *calculate_concs(paras)***

### **Input parameters:** 
`paras` : the class 'dict' object from **Set parameters** files (e.g., Start_SetParam_SA_example.py)

### **Output parameters:** 
`O2conc` : np.transpose([`O2conc1`,`O2conc2`]) \
`SO2conc` : np.transpose([`SO2conc1`,`SO2conc2`]) \
`H2Oconc` : np.transpose([`H2Oconc1`,`H2Oconc2`]) \
`paras` : updated class 'dict' object \
`export_file_folder` : the real export folder after conditional statement

### **Main parameters in the function:** 
**Notice:** the function will automatically determine the type of `flag_tube` ('4', '3', '2', '1' refer to the setup of the experiment)\
`flag_tube=4` : two tubes with different inner diameters and have Y piece, run the second tube with two flows simultaneously \
`flag_tube=3` : same as '4', but run the second tube with one flow after converting the mean concentrations \
`flag_tube=2` : two tubes with different inner diameters \
`flag_tube=1` : one tube

####  **flag_tube=1 & 2**
**Note that all the following flows and Q are needed in the input file but H2Oconc and T is optional.(See example input file under folder named Input_files)** \
`N2flow` : input N2 flow at different stages, unit sccm \
`O2flow` : input synthetic air flow at different stages, unit sccm \
`SO2flow` : input SO2 flow at different stages, unit sccm \
`H2Oflow` : input H2O flow at different stages, unit sccm \
`Q` : total flow in the tube at different stages, unit sccm (should be the same as sampleflow) \
`T` : temperature at different stages if recorded, unit K \
`H2O_concentration` : H2O concentration by measurement if recorded ('H2Oconc' at input file)

#### **flag_tube=3 & 4**
**Note that all the following flows and Q are needed in the input file but H2Oconc and T is optional.(See example input file under folder named Input_files)** \
`N2flow1`, `N2flow2` : input N2 flows for the first tube and Y piece at different stages, unit sccm \
`O2flow1`, `O2flow2` : input synthetic air flows for the first tube and Y piece at different stages, unit sccm \
`SO2flow` : input SO2 flow for first tube at different stages, unit sccm \
`H2Oflow1`, `H2Oflow2` : input H2O flows for the first tube and Y piece at different stages, unit sccm \
`Q1`, `Q2` : total flow in the first tube and the second tube (not Y piece!!!) at different stages, unit sccm (`Q2` should be the same as sampleflow) \
`T` : temperature at different stages if recorded, unit K \
`H2O_conc1`, `H2O_conc2` : H2O concentration for the first tube and Y piece by measurement if recorded ('H2Oconc1' and 'H2Oconc2' at input file)

#### **flag_tube=1 & 2 & 3 & 4**
`totalFlow1`, `totalFlow2` : total flow of calculating concentrations of gases in the first tube and the second tube at different stages, unit sccm \
`O2conc1`, `O2conc2` : calculated concentrations of O2 in the first tube and the second tube at different stages, unit cm-3 \
`SO2conc1`, `SO2conc2` : calculated concentrations of SO2 in the first tube and the second tube at different stages, unit cm-3 \
`H2Oconc1`, `H2Oconc2` : calculated concentrations of H2O in the first tube and the second tube at different stages, unit cm-3 \
`csH2O` : absorption cross section of water vapor (default = 7.22e-20), unit cm2 \
`qyH2O` : quantum yield (default = 1) \
`H2Oconc_free` : (only for `flag_tube=4`) H2O concentration from Y piece in the second tube, unit cm-3 \
`O2conc_free` : (only for `flag_tube=4`) O2 concentration from Y piece in the second tube, unit cm-3 \
`OHconc` : OH concentration by photolysis of H2O in the first tube at different stages, unit cm-3

## Run model (Run_flowtube.py)

### **Function *Run_flowtube(paras, export_file_folder, const_comp_conc, Init_comp_conc, num_stage)***
### **Input parameters:**
`paras` : the class 'dict' object from **Set parameters** files (e.g., Start_SetParam_SA_example.py) \
`export_file_folder` : the folder where the csv and txt files containing export infomration locate \
`const_comp_conc` : the constant concentrations for species from **Set parameters files** (e.g., Start_SetParam_SA_example.py) \
`Init_comp_conc` : the initial concentrations for species from **Set parameters files** (e.g., Start_SetParam_SA_example.py) \
`num_stage` : also from **Set parameters files** (e.g., Start_SetParam_SA_example.py), default is `num_stage= paras['OHconc']` but can be set by user as a integer number



### **Main parameters in the function**
`meanconc_s` : output table with headers containing steady state concentrations of species in `paras['plot_spec']` at different stages \
`c` : concentrations of species in `paras['plot_spec']` at all steps of all stages.