from molmass import Formula
import re
import pandas as pd
import numpy as np
import csv
def diff_coef(Mole, T, P, carrier_type,self):
    '''
    This function calculate diffusion coefficient according to Fuller's method (Fuller et al. 1965, 1966, 1969)
    Mole: molecular formula
    T: temperature in K
    P: pressure in pascal
    carrier_type: the type of carrier gas. The value should be N2, O2 or air
    '''
    if Mole in ['NA','M','MO2']:
        return [0,0,0]
    else:
        # if self.flag_mech == '1':
        #     data = pd.read_csv(self.tsv_file, sep='\t')
        #     data = data.dropna(subset=['Formula'])
        # else:
        #     data = pd.read_csv('SMILE_formula.csv', sep=';')
        # Detect delimiter
        with open(self.tsv_file, 'r', encoding='utf-8') as f:
            first_line = f.readline()  # Read first line
            dialect = csv.Sniffer().sniff(first_line)  # Detect delimiter

        # Read CSV with detected delimiter
        data = pd.read_csv(self.tsv_file, sep=dialect.delimiter)
         
        if carrier_type == "air":
            m_B = 29 #g mol-1
            V_B = 19.7 #no unit
        elif carrier_type == "O2":
            m_B = 32 #g mol-1
            V_B = 16.3 #no unit
        elif carrier_type == "N2":
            m_B = 28 #g mol-1
            V_B = 18.5  # no unit
        #calculate dimensionless diffusion volumes of A and B
        # atomic diffusion volume table. Fuller 1969
        Diff_vol = pd.DataFrame({'C': [15.9],
                                'H': [2.31],
                                'D': [2.31],
                                'O': [6.11],
                                'N': [4.54],
                                'F': [14.7],
                                'Cl': [21.0],
                                'Br': [21.9],
                                'I': [29.8],
                                'S': [22.9]})
        s = re.sub
        atom_count_mole = pd.DataFrame()
        # Look up formula from TSV; if not found, use species name as formula
        try:
            Mole = data.loc[data['Name'] == Mole, 'Formula'].values[0]
        except (IndexError, KeyError):
            pass  # keep original Mole (species name) as formula
        try:
            f_mole = s("[()',]", '', str(eval(s(',?(\d+)', r'*\1,', s('([A-Z][a-z]*)', r'("\1",),', Mole))))).split()
            num_O = f_mole.count('O')
            for c in set(f_mole): atom_count_mole[c] = [f_mole.count(c)]
        except:
            num_O = 0

        #calculate D for Mole
        V_A = np.array([0])
        for a in atom_count_mole.columns:
            V_A = V_A + atom_count_mole[a][0] * Diff_vol[a][0]
        #calculate molecular mass of the examined species
        m_A = Formula(Mole).mass
        m_AB = 2 / (1 / m_A + 1 / m_B)
        D_torr = 1.0868 * T ** 1.75 / (m_AB ** 0.5 * (V_A ** (1/3) + V_B ** (1/3)) ** 2) # torr cm2 s-1
        D = D_torr / (P / 133.3)

        return [float(D),m_A, num_O]
