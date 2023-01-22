# MARFORCE-Flowtube
This is a project to develop a flow tube chemistry module that include 2D space and chemistry processes

Shall we include the instruction for SA calibration experiments

# Table of contents
**[1. Documentation](#1-documentation)**

**[2. Installation](#2-installation)**

**[3. Running](#3-running)**

**[4. Inputs](#4-Inputs)**

* [4.1. Chemical Scheme file](#41-Chemical-Scheme-file)
* [4.2. SA calibration](#42-SA-calibration)
    + [4.2.1. flag_tube=1](#421-Inputs-for-flag_tube1)
    + [4.2.2. flag_tube=2](#422-Inputs-for-flag_tube2)
    + [4.2.3. flag_tube=3](#423-Inputs-for-flag_tube3)
    + [4.2.4. flag_tube=4](#424-Inputs-for-flag_tube4)
* [4.3. HOI calibration](#43-HOI-calibration)

**[5. Outputs](#5-outputs)**

**[6. Acknowledgements](#6-Acknowledgements)**

----

## 1. Documentation<a name="1-documentation"></a>
The README file you are now viewing serves as the MARFORCE (Marine Atmospheric paRticle FORmation and ChEmistry) manual, including how to run the model with correct inputs. 

The [article](NEED TO BE ADDED later) published in AMT explains the mechianisms of the flowtube model with corresponding schematics and simulation results. 

## 2. Installation<a name="installtion"></a>
1. Download the MARFORCE repository from https://github.com/momo-catcat/PANDA520-flowtube.
2. Create a environment containing Python and relevant libraries (). Anaconda is recommended to manage and install different libraries. 

## 3. Running<a name="running"></a>
1. After downloading the model package, go to *PANDA520-flowtube/PANDA520_flowtube/*
2. For model inputs, you should have: a .txt chemical reaction scheme file under folder *input_mechanism/* (e.g., 'SO2_SA.txt' given for SA calibration), a .csv file containing information of flows and possible temperature or concentrations under folder *Input_files/* (e.g., 'SA_cali_2021-09-10.csv') or any input folder set by yourself, a .py file to set the experimental information under the current folder (Start_SetParam_SA_example.py given as an example file for SA calibration). -- See next section **[4. Inputs](#4-Inputs)** for details. 
3. Once all the three files mentioned above set properly according to experiments, activate the environment containing all the packages, then run the .py file (e.g., Start_SetParam_SA_example.py) to get the model starting.
4. Finally, the model will show surface plots of all the steps for each experiment stage. The output files (one .csv and one.txt) will be saved under the folder *Export_files/* if no other folder is set for output. 

## 4. Inputs<a name="inputs"></a>
As mentioned above, there are three input files. This section will introduce how to set the files based on different mechanisms. 

### 4.1. Chemical Scheme file<a name="41-Chemical-Scheme-file"></a>
The chemical scheme file includes the reactions and their rate coefficients in the gas- and aqueous-phases.

Two example chemical scheme files are given under *PANDA520-flowtube/PANDA520_flowtube/input_mechanism*, named 'SO2_SA.txt' and 'HOI_cali_chem.txt' for SA and HOI calibration system. The chemical mechanistic information was taken from the Master Chemical Mechanism, MCM v3.3.1., via [website](http://mcm.york.ac.uk/). 

Markers are required to recognise different sections of the chemical scheme. The markers are for the MCM KPP format.

The expression for the rate coefficient can use Fortran type scientific notation or python type; acceptable math functions: EXP, exp, dsqrt, dlog, LOG, dabs, LOG10, numpy.exp, numpy.sqrt, numpy.log, numpy.abs, numpy.log10; rate coefficients may be functions of TEMP, RH, M, N2, O2 where TEMP is temperature (K), RH is relative humidity (0-1), M, N2 and O2 are the concentrations of third body, nitrogen and oxygen, respectively (# molecules/cc (air)). (Adopted from http://github.com/simonom/PyCHAM)

### 4.2. SA calibration<a name="42-SA-set-parameters"></a>
The chemical scheme file for sulfuric acid (SA) calibration is already given in the file named 'SO2_SA.txt', under *PANDA520-flowtube/PANDA520_flowtube/input_mechanism/*. 

**Notice:** if you use the given function file 'Calcu_by_flow_SA.py' under *PANDA520-flowtube/PANDA520_flowtube/Funcs/* for concentration calculation, the type of `flag_tube` will be determined automatically ('4', '3', '2', '1' refer to the setup of the experiment) according to the input inforamtion stating below. See the [article](NEED TO BE ADDED later) for schematics of different setups 
| Type of the experiment | Description of the flowtube|
|-------------|-------------|
|flag_tube=4 | Two tubes with different inner diameters and have Y piece, run the second tube with two flows simultaneously |
|flag_tube=3 | Same as '4', but run the second tube with one flow after converting the mean concentrations |
|flag_tube=2 | Two tubes with different inner diameters |
|flag_tube=1 | One tube |

#### **4.2.1. Inputs for flag_tube=1**<a name="421-Inputs-for-flag_tube1"></a>

**Flow variables .csv file:** an example is provided under folder *PANDA520-flowtube/PANDA520_flowtube/Input_files/*, called 'SA_cali_2021-09-10.csv'. It must include flow rate information but there are a few optional input variables. The headers show the variable name and the rest rows mean the number of the experiment stages. The model does not check whether the UV light is on or not, **thus just input all the stages with lights on or add your own codes while calculating the gas concentrations**.

| Input variables of the .csv file for flag_tube=1 | Description |
|-------------|-------------|
| N2flow | Input N2 flow at different stages (sccm)|
| O2flow | Input synthetic air flow at different stages (sccm) |
| SO2flow | Input SO2 flow at different stages (sccm)|
| H2Oflow | Input H2O flow at different stages (sccm)|
| Q | Total flow in the tube at different stages (sccm) |
| T *(optional)* | Temperature at different stages if recorded (K) |
| H2Oconc *(optional)*| H2O concentration by measurement if recorded (cm-3). If the measured H2O concentrations is more accurate than calculated ones by flows, the measured ones should be used for running the flowtube model. |


**Experimental information .py file:** an example is provided under folder *PANDA520-flowtube/PANDA520_flowtube/*, called 'Start_SetParam_SA_example.py'. A class 'dict' object is used for storing information.
| Input variables in the 'dict' object (paras) of the .py file for flag_tube=1 | Description|
|-------------|-------------|
| p | Pressure under which the experiment is conducted (Pa) |
| T | Temperature under which the experiment is conducted (K) |
| R1 |Inner radius for the first tube (cm), in this case flag_tube=1, thus the first tube is the only tube. |
| L1 | Length for the first tube (cm) |
| Itx | It product value for the specific UV lamp used in the experiment|
| Qx | The flow rate at which Itx is determined (lpm) |
| outflowLocation | Outflow (exhaust) tube located 'before' or 'after' injecting synthetic air, water vapor and SO2 |
| fullOrSimpleModel | 'simple' means Gormley & Kennedy approximation, while 'full' means flow model (much slower) |
| sampleflow |Inlet flow (lpm) of CIMs (chemical ionization mass spectrometer); it should be the same as total flow in the tube (Q). |
| SO2ratio | SO2 ratio of the gas bottle (ppm) |
| O2ratio | O2 ratio in synthetic air |
| Zgrid | Number of grids in direction of tube length  |
| Rgrid | Number of grids in direction of radius |
| dt | Differential time interval (usually 1e-4, but try small number if the values are too large out of range and showing 'Nan' during calculation) |
| model_mode | Use 'normal' if you don't know what this is for. 'kinetic' mode refers to running the model without chemistry module to test the kinetic core. |
| Diff_setname | Diffusion for the species that you want to define by yourself, otherwise it will be calculated automatically based on the elements it contains |
| Diff_set | Add the value according to the Diff_setname |
| sch_name | Chemical scheme file name stored in the *PANDA520-flowtube/PANDA520_flowtube/input_mechanism/* folder|
| const_comp | Species you think they should have constant concentrations in the whole tube (for SA calibration, const_comp=['SO2','O2','H2O']) |
| Init_comp | Species you think they should have initial concentration in the first grid of tube (for SA calibration, Init_comp=['OH','HO2'])  |
| key_spe_for_plot | Key species as criterion to stop the loop (for SA calibration, key_spe_for_plot='H2SO4') |
| plot_spec | Species that you want to plot |
| file_name | The name of the flow variables .csv file. |
| input_file_folder *(optional)* | The folder where the flow variables .csv file locates; default folder is *PANDA520-flowtube/PANDA520_flowtube/Input_files/* if this parameter is missing  |
| export_file_folder *(optional)* | The folder where the .csv and .txt files containing export infomration (such as H2SO4 concentration) locate; default folder is *PANDA520-flowtube/PANDA520_flowtube/Export_files* if this paramter is missing |
| flag_tube *(optional)* | The type of experiment setup, it will determined automatically if 'Calcu_by_flow_SA.py' is used for concentration calculation|

For the other input variables of the .py file stated below, **you don't need to change if** you use the example file 'Start_SetParam_SA_example.py' with the 'Calcu_by_flow_SA.py' as the function to calculate gas concentrations. **You could also make your own script for calculating concentrations based on flows, and in this case you also have to determine the flag_tube by yourself. We RECOMMEND you to use the example function file 'Calcu_by_flow_SA.py'**  
| Other input variables of the .py file | Description|
|-------------|-------------|
| Init_comp_conc | Initial concentrations for species that you already set in the paras (in the case of SA calibration, 'OHconc' (concentration of OH radical calculated from Itx, Qx, H2O concentration and etc.)  is used for both 'OH' and 'HO2') |
| const_comp_conc | Constant concentrations for species you already set |
| num_stage | Default is num_stage= paras['OHconc'], meaning the number of stages based on OHconc since we want to calculate the stages with light on, but **this can also be set by user as a integer number N (in this case, it means the first N stages are calculated)** |


#### **4.2.2. Inputs for flag_tube=2**<a name="422-Inputs-for-flag_tube2"></a>
**Flow variables .csv file** is set exactly the same as the the situation of **flag_tube=1**, see details above. \
For  **experimental information .py file**, all the variables are the same except here R2 and L2 are additionally needed.
| Additional input variables in the 'dict' object (paras) of the .py file for flag_tube=2 | Description|
|-------------|-------------|
| R2 | Inner diameter for the second tube (cm) if there is any, it could also be set to 0 for flag_tube = 1 |
| L2 | Length for the second tube (cm) |

#### **4.2.3. Inputs for flag_tube=3**<a name="423-Inputs-for-flag_tube3"></a>
**Flow variables .csv file:** since a Y piece is added into the setup, variables for flows coming from Y piece is needed to input.

| Input variables of the .csv file for flag_tube=3 | Description |
|-------------|-------------|
| N2flow1 | Input N2 flow for the first tube at different stages (sccm)|
| N2flow2 | Input N2 flow for the Y piece at different stages (sccm)|
| O2flow1 | Input synthetic air flow for the first tube at different stages (sccm) |
| O2flow2 | Input synthetic air flow for the Y piece at different stages (sccm) |
| SO2flow | Input SO2 flow for the first tube at different stages (sccm)|
| H2Oflow1 | Input H2O flow for the first tube at different stages (sccm)|
| H2Oflow2 | Input H2O flow for the Y piece at different stages (sccm)|
| Q1 | Total flow in the first tube at different stages (sccm) |
| Q2 | Total flow in the second tube (NOT Y piece!!!) at different stages (sccm) |
| T *(optional)* | Temperature at different stages if recorded (K) |
| H2Oconc1 *(optional)*| H2O concentration by measurement for the first tube if recorded (cm-3). If the measured H2O concentrations is more accurate than calculated ones by flows, the measured ones should be used for running the flowtube model. |
| H2Oconc2 *(optional)*| H2O concentration by measurement for the Y piece if recorded (cm-3).|

**Experimental information .py file:** variables are the same as those for flag_tube=2 if 'Calcu_by_flow_SA.py' is used for calculating concentrations, otherwise you should make your own script for calculating the concentrations in both the first and second tubes. 


#### **4.2.4. Inputs for flag_tube=4**<a name="424-Inputs-for-flag_tube4"></a>
Every input variables are the same as those for flag_tube=3, except that for **Experimental information .py file**, you have to set flag_tube=4 by yourself while using the example script 'Calcu_by_flow_SA.py'. 

| Additional input variables in the 'dict' object (paras) of the .py file for flag_tube=4 | Description|
|-------------|-------------|
| flag_tube | In this case ('4'), this variable is required to write by yourself, but if you want to make your own script, this is not required. |

### 4.3. HOI calibration<a name="43-HOI-set-parameters"></a>

## 5. Outputs<a name="5-outputs"></a>

## 6. Acknowledgements<a name="6-Acknowledgements"></a>



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

