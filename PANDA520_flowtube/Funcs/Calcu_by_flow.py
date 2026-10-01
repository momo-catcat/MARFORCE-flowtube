import numpy as np
import pandas as pd
import os

from pandas import Index
from Funcs.Vapour_calc import H2O_conc
from Funcs.convert_tsv2xml import tsv_to_custom_xml


def calculate_concentrations(modelparams):
    """
    Calculate gas concentrations and other parameters for the flow tube experiment.

    Parameters:
        modelparams (SimulationParams): Object containing experimental parameters.

    Returns:
        SimulationParams: Updated modelparams with calculated concentrations and paths.
    """
    # Define paths
    dirpath = os.path.dirname(__file__) #if '__file__' in globals() else os.getcwd()  # Get current file path
    input_file_folder = getattr(modelparams, 'input_file_folder', os.path.normpath((os.path.join(dirpath, '../Input_files/'))))
    export_file_folder = getattr(modelparams, 'export_file_folder', os.path.normpath((os.path.join(dirpath, '../Export_files/'))))
    input_mechanism_folder = os.path.normpath(os.path.join(getattr(modelparams, 'input_mechanism_folder', os.path.join(dirpath, '../input_mechanism/')), modelparams.folder_mechaism))

    subfolder_name = modelparams.file_name[:-4]
    export_subfolder_path = os.path.join(export_file_folder, subfolder_name)

    # Make sure the directory exists
    os.makedirs(export_subfolder_path, exist_ok=True)
    
    # check if it is a point source of OH or not 
    if hasattr(modelparams, 'Rl') and (modelparams.Rl > 0):
        modelparams.OHsource = 'Continuous'
    else:
        modelparams.OHsource = 'point'

    # Load input data
    data = pd.read_csv(os.path.join(input_file_folder, modelparams.file_name))

    # Identify flow and concentration names
    flowname = [i.replace('flow', '_') if i[-1].isdigit() else i.replace('flow', '') for i in data.columns.to_list() if 'flow' in i]
    flowname_ori = [i for i in data.columns.to_list() if 'flow' in i]

    concname = [i.replace('conc', '_') if i[-1].isdigit() else i.replace('conc', '') for i in data.columns.to_list() if 'conc' in i]
    concname_ori = [i for i in data.columns.to_list() if 'conc' in i]

    # Extract parameters
    paras_ratios = [key.replace('ratio', '') for key in modelparams.__dict__ if ('ratio' in key) ]
    paras_ratios.remove('O2')
    modelparams.runtime = data['time'].values if 'time' in data.columns else np.arange(len(data)) * 0
    outflow_location = modelparams.outflowLocation
    sample_flow = modelparams.sampleflow
    pressure = modelparams.p
    temperature = modelparams.TEMP
    Itx = modelparams.Itx
    Qx = modelparams.Qx

    inputs: Index = data.columns.tolist()
    # Determine flow tube type
    basic_flow_columns = ['H2Oflow', 'N2flow']
    flag_tube = None
    modelparams.actual_acc = modelparams.timesteps * modelparams.inter_acc/10000
    if all(col in inputs for col in basic_flow_columns):
        # Single-tube setup (Tube 1 or 2)
        flows = {col: data[orig_col] for col, orig_col in zip(flowname, flowname_ori)} 
        concs = {col: data[orig_col] for col, orig_col in zip(concname, concname_ori)} 
        concs = {key: value.dropna() for key, value in concs.items() if not value.dropna().empty}

        flows['Q'] = data['Q']  # Add 'Q' to the dictionary separately
        sum_flow = sum(flows[i] for i in flowname)
        ion_trap_intensity = Itx * Qx / (flows['Q'] / 1000)
        idx = flows['H2O'] < 1000

        # Check if Tube 1 or 2
        if hasattr(modelparams, 'L2'):
            flag_tube = '2'
        else:
            flag_tube = '1'
            modelparams.L2 = np.float64(0)
            modelparams.R2 = modelparams.R1
        ## get all the transponsed concentration 
        conc_groups = concs.keys()
        transposed_dict = {}
        for group, columns in concs.items():
            transposed_dict[group] = np.transpose([columns, columns])

    else:
        # Dual-tube setup (Tube 3)
        flows = {col: data[orig_col] for col, orig_col in zip(flowname, flowname_ori)} 
        concs = {col: data[orig_col] for col, orig_col in zip(concname, concname_ori)} 
        concs = {key: value.dropna() for key, value in concs.items() if not value.dropna().empty}
        flows['Q1'] = data['Q1']  # Add 'Q' to the dictionary separately
        flows['Q2'] = data['Q2']  # Add 'Q' to the dictionary separately

        cols1 = [i for i in flowname if '2' not in i] 
        sum_flow = sum(flows[i] for i in cols1)
        ion_trap_intensity = Itx * Qx / (flows['Q1'] / 1000)
        flag_tube = '3'
        idx = flows['H2O_1'] < 1000

        ## get all the transponsed concentration 
        conc_groups = {}

        for key in concs.keys():
            # Split the key to get the prefix (e.g., 'H2O', 'X', 'Y')
            prefix = key.split('_')[0]
            
            # Add the column to the appropriate group
            if prefix not in conc_groups:
                conc_groups[prefix] = []
            conc_groups[prefix].append(concs[key])

        # Transpose the selected columns for each group (H2O, X, Y, etc.)
        transposed_dict = {}
        for group, columns in conc_groups.items():
            transposed_dict[group] = np.transpose([col for col in columns])

  # Handle temperature data
    temperature = data['T'] if 'T' in data.columns else np.ones(len(data)) * temperature

    kB = 1.3806488e-23  # Boltzmann constant

    # Calculate concentrations based on flow tube type
    if flag_tube == '3':
        total_flow_1 = sum_flow if outflow_location == 'after' else flows['Q1']
        total_flow_2 = sample_flow * 1e3
        H2O_conc_1 = [flows['H2O_1'][i] / total_flow_1[i] * H2O_conc(temperature[i], 1).SatP[0] / kB / temperature[i] / 1e6 for i in range(len(temperature))]
        O2_conc_1 = [flows['O2_1'][i] * modelparams.O2ratio / total_flow_1[i] * pressure / kB / temperature[i] / 1e6 for i in range(len(temperature))]
        modelparams.Q1, modelparams.Q2 = flows['Q1'], flows['Q2']
    else:
        total_flow_1 = sum_flow if outflow_location == 'after' else flows['Q']
        H2O_conc_1 = [flows['H2O'][i] / total_flow_1[i] * H2O_conc(temperature[i], 1).SatP[0] / kB / temperature[i] / 1e6 for i in range(len(temperature))]
        O2_conc_1 = [flows['O2'][i] * modelparams.O2ratio / total_flow_1[i] * pressure / kB / temperature[i] / 1e6 for i in range(len(temperature))]
        modelparams.Q1, modelparams.Q2 = flows['Q'], flows['Q']

    ### calculate the concentration of other species expect O2
    for sp in paras_ratios:
        try:
            flows[sp+'_conc_1'] = [ flows[sp][i] * getattr(modelparams, sp+'ratio') / total_flow_1[i] * pressure / kB / temperature[i] / 1e6 for i in range(len(temperature))]
        except:
            pass

    if flag_tube in ['1', '2']:
        H2O_conc_2, O2_conc_2 = H2O_conc_1, O2_conc_1
        for sp in paras_ratios:
            try:
                flows[sp+'_conc_2'] =flows[sp+'_conc_1']
            except:
                pass
    else:
        total_flow_2 = sample_flow * 1e3
        H2O_conc_2 = [(H2O_conc_1[i] * total_flow_1[i] + flows['H2O_2'][i] * H2O_conc(temperature[i], 1).SatP[0] / kB / temperature[i] / 1e6) / total_flow_2 for i in range(len(temperature))]
        O2_conc_2 = [(O2_conc_1[i] * total_flow_1[i] + flows['O2_2'][i] * modelparams.O2ratio * pressure / kB / temperature[i] / 1e6) / total_flow_2 for i in range(len(temperature))]
        for sp in paras_ratios:
            try:
                flows[sp+'_conc_2'] = [flows[sp+'_conc_1'][i] * total_flow_1[i] / total_flow_2 for i in range(len(flows[sp+'_conc_1']))]
            except:
                pass

    # Calculate OH concentration
    OH_conc = ion_trap_intensity * 7.22e-20 * 1 * np.array(H2O_conc_1)
    
    # determine the H2O finally 
    if 'H2O' in conc_groups:
        H2Oconc = transposed_dict['H2O']
    else:
        H2Oconc = np.transpose([H2O_conc_1, H2O_conc_2])
    # use the calucated H2O when H2O flow is below 1000
    H2Oconc[idx] = np.transpose([np.array(H2O_conc_1)[idx], np.array(H2O_conc_2)[idx]])

    if 'flag_tube' in modelparams.__dict__:
        flag_tube1 = modelparams.flag_tube
        if flag_tube1 == '4' and flag_tube == '3':
            flag_tube = flag_tube1

    if flag_tube == '4':
        const_comp_free = ['H2O', 'O2']
        H2Oconc_free = [flows['H2O_2'][i] / (total_flow_2-total_flow_1[i]) * H2O_conc(temperature[i], 1).SatP[0] / kB / temperature[i] / 1e6 for i in range(len(temperature))]
        O2conc_free = [flows['O2_2'][i] / (total_flow_2-total_flow_1[i]) *modelparams.O2ratio * pressure / kB / temperature[i] / 1e6  for i in
                   range(len(temperature))]
        const_comp_conc_free = [H2Oconc_free, O2conc_free]
    else:
        const_comp_free = []
        const_comp_conc_free = [0]
    
    # Update modelparams
    modelparams.O2conc = np.transpose([O2_conc_1, O2_conc_2])
    modelparams.H2Oconc = H2Oconc
    modelparams.OHconc = OH_conc
    modelparams.HO2conc = OH_conc

    modelparams.flag_tube = flag_tube
    modelparams.sch_name = os.path.join(input_mechanism_folder, modelparams.sch_name)
    modelparams.tsv_file = os.path.join(input_mechanism_folder, modelparams.tsv_file)

    modelparams.export_file_folder = export_subfolder_path +'/'
    modelparams.num_sp_out_O2 = len(paras_ratios)
    modelparams.const_comp_free = const_comp_free
    modelparams.const_comp_conc_free = const_comp_conc_free

    for sp in paras_ratios:
        try:
            setattr(modelparams, sp + 'conc', np.transpose([flows[sp + '_conc_1'], flows[sp + '_conc_2']]))
        except KeyError:
            pass  # Handles cases where a species concentration is missing
    # remove water from the transposed_dict
    transposed_dict.pop('H2O', None)
    # print('transposed_dict =', transposed_dict)
    for group, transposed in transposed_dict.items():
        if group in modelparams.const_comp:
            transposed = np.asarray(transposed).flatten()
            transposed2 = transposed / sample_flow * (modelparams.Q1/1e3)
            combined = np.column_stack([transposed, transposed2])
            setattr(modelparams, group + 'conc', combined)
        else:
            setattr(modelparams, group + 'conc', np.asarray(transposed).flatten())

    xml_file = "chemical_species_custom.xml"  # Desired output XML file path

    if modelparams.flag_mech == '1':
        # Call the function
        tsv_to_custom_xml(modelparams.tsv_file, input_mechanism_folder, xml_file)
    else: # give a default xml file and this is not going to work for other mechanism
        pass
    
    modelparams.xml_name = os.path.join(input_mechanism_folder, 'chemical_species_custom.xml')
    # add the constant compounds concentration
    const_comp_conc = np.transpose([getattr(modelparams, k)  for k in [i + 'conc' for i in modelparams.const_comp]])
    modelparams.const_comp_conc = const_comp_conc
    Init_comp_conc = [getattr(modelparams, k)  for k in [i + 'conc' for i in modelparams.Init_comp]]
    try:
        modelparams.Init_comp_conc = np.column_stack(Init_comp_conc)
    except ValueError:
        pass
    modelparams.funcs_folder = dirpath
    modelparams.photo_path = os.path.join(dirpath, '../photofiles/MCMv3.2/')

    return modelparams

