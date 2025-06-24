# Functions and utils for turbo-mini-space-cowboy notebook
# Basile Wicky -- 240116
# v8- updated concentration measurements to not fail if sample is too concentrated to normalize -- srgerb 240822

# Libraries
import glob
import os
import json
from datetime import date; today=date.today().strftime('%Y%m%d')[2:]
import numpy as np
import re
from pycorn import pc_uni6
from Bio import SeqIO, Seq
from scipy import stats, signal
from scipy.cluster import hierarchy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import getpass

##################
# ECHO stuff
##################

def gen_echo_transfers(
    eblock_plates,
    vectors=['LM0627'],
    vector_stock_conc=[100],
    vector_amount_fmol=4,
    path_to_vectors='/software/lab/johnbercow/entry_vectors/'
):
    '''
    Generate ECHO transfer file. 
    
    NB: Currently implemented for cloning every specified eBlock into every sepecifed vector.
    
    eblocks_plates: dictionary containing path to .xlsx files from IDT eBlocks order(s).
    
    vectors: list of GG entry vector IDs. Each eBlock will be cloned into each vector.
        If a 'Vector' column is included in the spreadsheet, this information will be parsed instead.
        If specifying vectors in the spreadsheet, then all vectors used (and their concentrations) need to be specified. 
    
    vector_stock_conc: list of concentrations for the GG entry vector stocks, in ng/uL.
    
    vector_amount_fmol: amount of vector to use for the reaction, in fmol (default: 4 fmol). 
        Use less if you cannot get a positive H2O volume value in the master mix.
    
    path_to_vectors: path to directory containing entry vectors FASTA files. 
    '''

    w96 = [r + str(c) for r in 'ABCDEFGH' for c in range(1, 13)]    
    w384 = [r + str(c) for r in 'ABCDEFGHIJKLMNOP' for c in range(1, 25)]
    
    if len(vectors) != len(vector_stock_conc):
        print('Specify the concentration of each vector stock!\nOnly entries for vectors with specified concentrations will be generated.') 
    
    # Hard-coded defintions.
    tot_rxn_vol = 1 # uL
    eblock_conc = 4 # ng/uL, as delivered by IDT in ECHO-qualified plates
    insert_to_vector = 2 # insert:vector molar ratio
    
    # Sources.
    echo_df = pd.DataFrame()
    for k, v in eblock_plates.items():
        tmp =  pd.read_excel(v).rename(
            columns={'Sequence':'eblock', 'Well Position':'Source Well'}
        ).dropna()
        tmp['Source Plate Name'] = k
        echo_df = pd.concat([echo_df, tmp])
    
    echo_df.reset_index(drop=True, inplace=True)
    
    
    # Destinations (assumes 348w (eBlocks) --> 96w (cultures)).
    split_idx = [echo_df.index.values[x:x+96] for x in range(0, len(echo_df), 96)]
    for i, si in enumerate(split_idx):
        echo_df.loc[si, 'Destination Plate Name'] = i + 1
        echo_df.loc[si, 'Destination Well'] = w96[:len(si)]
    
    
    # eBlock transfer volume.
    echo_df['eblock_length'] = echo_df['eblock'].str.len()
    eblock_bp = echo_df['eblock_length'].median()   
    eblock_vol = ((vector_amount_fmol * insert_to_vector) / (((eblock_conc * 1e-9) / ((eblock_bp * 617.96) + 36.04)) * 1e15))
    echo_df['Transfer Volume'] = eblock_vol * 1e3 # transfer volumes are specified in nL


    tmps = {}
    if 'Vector' in echo_df.columns: # clone into row-specific vectors if they are included in the spreadsheet.
        for vector, sub_echo_df in echo_df.groupby('Vector'):
            tmps[vector] = sub_echo_df
    
    else: # clone each eBlock into each vector specified with the function if the spreadsheet does not contain a 'Vector' column.
        for i, (vector, vector_conc) in enumerate(list(zip(vectors, vector_stock_conc))):
            tmp = echo_df.copy()
            tmp['Vector'] = vector
            n_destp = len(echo_df['Destination Plate Name'].unique())
            tmp.loc[:, 'Destination Plate Name'] = tmp['Destination Plate Name'] + (i * n_destp)
            tmps[vector] = tmp
        
    mm_dfs = []
    echo_dfs = []
    mm_wells = 0
    for i, (vector, vector_conc) in enumerate(list(zip(vectors, vector_stock_conc))):
        tmp = tmps[vector]

        # Master mix.
        mm_df = pd.DataFrame()

        # Per reaction calculations.
        vector_bp = len([str(r.seq) for r in SeqIO.parse(glob.glob(f'{path_to_vectors}{vector}*.fa')[0], 'fasta')][0])
        vector_vol = vector_amount_fmol / (((vector_conc * 1e-9) / ((vector_bp * 617.96) + 36.04)) * 1e15)
        mm_df.loc['T4 buffer', 'per_rxn_uL'] = tot_rxn_vol / 10
        mm_df.loc[f'{vector} @ {vector_conc} ng/uL', 'per_rxn_uL'] = vector_vol
        mm_df.loc['BsaI_HFv2 (NEB#R3733L)', 'per_rxn_uL'] = 0.06 # 1.2 units, from commercial stock @ 20 U/uL
        mm_df.loc['T4 ligase (NEB#M0202L)', 'per_rxn_uL'] = 0.1 # 40 units, from commercial stock @ 400 U/uL
        int_vol = mm_df.per_rxn_uL.sum()
        mm_df.loc['H2O', 'per_rxn_uL'] = tot_rxn_vol - eblock_vol - int_vol
        mm_df = mm_df.reindex([
            'H2O', 
            'T4 buffer', 
            f'{vector} @ {vector_conc} ng/uL', 
            'BsaI_HFv2 (NEB#R3733L)', 
            'T4 ligase (NEB#M0202L)'], 
        )

        # Volume correction to account for glycerol content in master mix.
        perc_glycerol = 50 * (mm_df.loc['BsaI_HFv2 (NEB#R3733L)', 'per_rxn_uL'] + mm_df.loc['T4 ligase (NEB#M0202L)', 'per_rxn_uL']) / mm_df.per_rxn_uL.sum()
        mm_corr = 0.49 * perc_glycerol - 0.008 * perc_glycerol**2 + 1 # correction based on calibration curve


        # Totals.
        n_rxns = len(tmp)
        mm_vol_per_well = tot_rxn_vol - eblock_vol
        mm_used_tot = mm_vol_per_well * n_rxns
        mm_needed_tot = (20 + mm_corr) * np.ceil(mm_used_tot / (65 - 20 - mm_corr)) + mm_used_tot
        mm_df.loc[:, 'tot_uL'] = mm_df.per_rxn_uL * n_rxns * mm_needed_tot / mm_used_tot 
        
        tmp2 = tmp[['Destination Plate Name', 'Destination Well', 'Vector']].reset_index(drop=True)
        tmp2['Source Plate Name'] = list(eblock_plates.keys())[-1]
        tmp2['Transfer Volume'] = mm_vol_per_well * 1e3 # transfer volumes are specified in nL

        # Assign MM wells.
        n_mm_srcw = np.ceil(mm_df.tot_uL.sum() / 65).astype(int)
        vol_mm_per_srcw = mm_df.tot_uL.sum() / n_mm_srcw
        n_destw = np.ceil(len(tmp) / n_mm_srcw).astype(int)
        for j in range(n_mm_srcw):
            tmp2.loc[j*n_destw:(j+1)*n_destw, 'Source Well'] = w384[::-1][mm_wells + j]
        
        mm_wells += n_mm_srcw
        
        mm_df.loc[f'Add [{vector}] MM to', 'per_rxn_uL'] = ', '.join(tmp2['Source Well'].unique())
        mm_df.loc[f'Add [{vector}] MM to', 'tot_uL'] = vol_mm_per_srcw
        mm_df.loc['', :] = ''
        mm_dfs.append(mm_df)
                
        echo_dfs.append(pd.concat([tmp, tmp2]))

        
    echo_df = pd.concat(echo_dfs).dropna()
    echo_df.loc[:, 'Destination Plate Name'] = echo_df['Destination Plate Name'].apply(lambda x: int(x))
    mm_df = pd.concat(mm_dfs)
    
    echo_df.to_csv(f'{today}_ECHO_cloning.csv', index=False)

    return echo_df, mm_df


#########################
# Cloning stuff
#########################


def clone(
    r,
    path_to_vectors='/software/lab/johnbercow/entry_vectors/'
):
    '''
    Generate plasmid map.
    
    path_to_vectors: path to directory containing entry vectors FASTA files. 
    '''

    # Extract entry vector sequences.
    gg_vectors = glob.glob(path_to_vectors + '*.fa')
    entry_vectors = {}
    for v in gg_vectors:
        records = SeqIO.parse(v, 'fasta')
        for record in records:
            entry_vectors[record.id.split('_')[0]] = str(record.seq.lower())

    # FW, RV, N_spacer, N_sticky -- for generating the plasmid maps
    cuts = {
        'BsaI':['ggtctc', 'gagacc', 1, 4],
        'SapI':['gctcttc', 'gaagagc', 1, 3]
    }
  
    vector = r['Vector']
    enzyme = (lambda x: 'SapI' if 'MA' in x else 'BsaI')(vector)

    # Enzyme-specific cut characteristics.
    fw, rv, n_spacer, n_sticky = cuts[enzyme]

    # GG cloning.
    entry_vector = entry_vectors[vector]
    vector_5prime, vector_3prime = entry_vector.find(fw) \
                                    + len(fw) \
                                    + n_spacer, entry_vector.find(rv) \
                                    - n_spacer
    
    insert = r['eblock'].lower()
    
    # Find eBlock section with cut sites facing in the correct directions.
    # (Necessary for cases where cut sites are accidentlly also present in the padding regions.)
    fw_locations = np.array([x.span()[0] for x in re.finditer(fw, insert)])
    rv_locations = np.array([x.span()[0] for x in re.finditer(rv, insert)])

    fw_idx = []
    rv_idx = []
    delta_bp = []
    for i, f in enumerate(fw_locations):
        for j, r in enumerate(rv_locations):
            delta_bp.append(r - f)
            fw_idx.append(i)
            rv_idx.append(j)
        
    delta_bp = np.array(delta_bp)
    correct_idx = np.argwhere(delta_bp==delta_bp[delta_bp>0].min())[0][0]

    insert_5prime = fw_locations[fw_idx[correct_idx]] \
                    + len(fw) \
                    + n_spacer \
                    + n_sticky

    insert_3prime = rv_locations[rv_idx[correct_idx]] \
                    - n_spacer \
                    - n_sticky
                        
    assembled_plasmid = entry_vector[:vector_3prime] \
                            + insert[insert_5prime:insert_3prime] \
                            + entry_vector[vector_5prime:]

    return assembled_plasmid.lower()


def coding_seq(
    row,
):
    '''
    Identify the full ORF that contains the insert.
    '''
    
    eblock = row['eblock'].lower() # the full eBlock sequence as ordered.

    # FW, RV, N_spacer, N_sticky -- for identifying the insert sequence.
    cuts = {
        'BsaI':['ggtctc', 'gagacc', 1, 4],
        'SapI':['gctcttc', 'gaagagc', 1, 3]
    }
  
    vector = row['Vector']
    enzyme = (lambda x: 'SapI' if 'MA' in x else 'BsaI')(vector)

    # Enzyme-specific cut characteristics.
    fw, rv, n_spacer, n_sticky = cuts[enzyme]
        
    # Find eBlock section with cut sites facing in the correct directions.
    # (Necessary for cases where cut sites are accidentlly also present in the padding regions.)
    fw_locations = np.array([x.span()[0] for x in re.finditer(fw, eblock)])
    rv_locations = np.array([x.span()[0] for x in re.finditer(rv, eblock)])

    fw_idx = []
    rv_idx = []
    delta_bp = []
    for i, f in enumerate(fw_locations):
        for j, r in enumerate(rv_locations):
            delta_bp.append(r - f)
            fw_idx.append(i)
            rv_idx.append(j)
        
    delta_bp = np.array(delta_bp)
    correct_idx = np.argwhere(delta_bp==delta_bp[delta_bp>0].min())[0][0]

    insert_5prime = fw_locations[fw_idx[correct_idx]] \
                    + len(fw) \
                    + n_spacer \
                    + n_sticky

    insert_3prime = rv_locations[rv_idx[correct_idx]] \
                    - n_spacer \
                    - n_sticky
                        
    insert = eblock[insert_5prime:insert_3prime] # the insert sequence without cloning adapters.

    plasmid_seq = row['plasmid_seq']

    # Identify the ORF that contains the insert.
    # Search for the shortest START-STOP span that contains the insert sequence.
    starts = np.array([s.start() for s in re.finditer('atg', plasmid_seq)])
    ends = np.array(sorted([e.end() for e in re.finditer('tag', plasmid_seq)] 
                           + [e.end() for e in re.finditer('taa', plasmid_seq)] 
                           + [e.end() for e in re.finditer('tga', plasmid_seq)]))
    current_stop = 0

    for s in starts:
        inframe_stops = ends[np.logical_and(ends>s, (ends-s)%3==0)]

        if len(inframe_stops) > 0:
            if s > current_stop:
                current_stop = inframe_stops[0]
                coding_seq = plasmid_seq[s:current_stop]
                
                if insert in coding_seq:
                    return coding_seq
    
    

############################
# Protparam stuff
############################

def e280(
    seq, 
    cystine=False,   
):
    '''
    Calculate the extinction coefficient at 280 nm (in M^-1 cm^-1).
    '''

    table_e280 = {
            'W':5500,
            'Y':1490,
            'S-S':125
            }

    e280 = 0
    for aa in 'WY':
        n_aa = seq.count(aa)
        e280 += n_aa * table_e280[aa]

    if cystine == True:
        e280 += (seq.count('C') / 2 ) * table_e280['S-S']
    
    return e280


def mw(
    seq,
):
    '''
    Calculate the molecular mass (in Da).
    '''

    table_avg_mass = {
            'A':71.0788,
            'R':156.1875,
            'N':114.1038,
            'D':115.0886,
            'C':103.1388,
            'E':129.1155,
            'Q':128.1307,
            'G':57.0519,
            'H':137.1411,
            'I':113.1594,
            'L':113.1594,
            'K':128.1741,
            'M':131.1926,
            'F':147.1766,
            'P':97.1167,
            'S':87.0782,
            'T':101.1051,
            'W':186.2132,
            'Y':163.1760,
            'V':99.1326,
            'H2O':18.01524,
            '*':0.0
            }

    mw = table_avg_mass['H2O']
    for aa in seq:
        mw += table_avg_mass[aa]

    return mw


def aliphatic_idx(
    seq,
):
    '''
    Calculate the aliphatic index.
    '''    
    
    percA = 100 * seq.count('A') / len(seq)
    percV = 100 * seq.count('V') / len(seq) 
    percI = 100 * seq.count('I') / len(seq)
    percL = 100 * seq.count('L') / len(seq)
    
    aliphatic_idx = percA + 2.9 * percV + 3.9 * (percI + percL)

    return aliphatic_idx


def charge_pI(
    seq, 
    pH_value=7.0,
):
    '''
    Calculate the isolectric point and the charge at the specified pH value (or 7.0 if unspecified).
    '''

    table_pKa = {
            'Nterm':8.23,
            'Cterm':3.55,
            'D':3.86,
            'E':4.34,
            'H':6.45,
            'C':8.49,
            'Y':9.76,
            'K':10.34,
            'R':13.9
            }

    aa_count = {aa:0 for aa in 'DEHCYKR'}
    for aa in seq:
        if aa in 'DEHCYKR':
            aa_count[aa] += 1

    charge_pH = {round(pH, 2):0 for pH in np.arange(0, 14.01, 0.01)}
    
    for pH in charge_pH.keys():
        charge = 0
        
        for aa in 'KRH':
            charge += aa_count[aa] / (1 + 10**(pH - table_pKa[aa]))
        
        for aa in 'DECY':
            charge -= aa_count[aa] / (1 + 1/(10**(pH - table_pKa[aa])))
         
        charge += 1 / (1 + 10**(pH - table_pKa['Nterm']))
        charge -= 1 / (1 + 1 / (10**(pH - table_pKa['Cterm'])))

        charge_pH[pH] = charge
   
    charge = charge_pH[pH_value]
    pI = min(charge_pH, key=lambda x: abs(charge_pH[x]))

    return charge, pI

def e205(
    seq, 
    cystine=False,
): 
    '''
    Calculate the extinction coefficient at 205 nm (in M^-1 cm^-1).
    '''

    table_e205 = {
            'W':20400,
            'F':8600,
            'Y':6080,
            'H':5200,
            'M':1830,
            'R':1350,
            'C':690,
            'N':400,
            'Q':400,
            'S-S':820,
            'peptide_bb':2780
            }

    # Values from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3690723/
    # NB1 error on peptide_bb is 2780 +/- 168
    # NB2 the S-S value is for the bond only, the individual contributions from each Cys (690) still need to be counted

    e205 = (len(seq) - 1) * table_e205['peptide_bb']
    
    for aa in 'WFYHMRCNQ':
        n_aa = seq.count(aa)
        e205 += n_aa * table_e205[aa]

    if cystine == True:
        e205 += (seq.count('C') / 2 ) * table_e205['S-S']

    return e205

def protparam(echo_df, oligomers, sec_buffer_pH):
    '''
    Add protoparam informations to dataframe.
    
    oligomers: Boolean variable indicating whether or not the design is expected to be an oligomer. 
        If True, a n_chains column needs to be specified in the dataframe.
    
    sec_buffer_pH: the pH of the buffer used for SEC. Influences charge calculations.
    '''

    # Obtain expression product from ORF.
    echo_df['exp_prod'] = echo_df['ORF'].apply(lambda x: str(Seq.Seq(x).translate()))

    # Compute protein paramters.
    echo_df['nAA'] = echo_df['exp_prod'].str.len()
    echo_df['protomer_MW'] = echo_df['exp_prod'].map(mw)
    echo_df['aliphatic_idx'] = echo_df['exp_prod'].map(aliphatic_idx)
    echo_df['e280'] = echo_df['exp_prod'].map(e280)
    echo_df['OD280'] = echo_df['e280'] / echo_df['protomer_MW']
    echo_df['e205'] = echo_df['exp_prod'].map(e205)
    echo_df['OD205'] = echo_df['e205'] / echo_df['protomer_MW']
    echo_df['nC'] = echo_df['exp_prod'].str.count('C')
    echo_df[f'charge@{sec_buffer_pH}'], echo_df['pI'] = zip(
        *echo_df['exp_prod'].apply(lambda x: charge_pI(x, sec_buffer_pH))
    )

    if oligomers:
        echo_df['oligomer_MW'] = echo_df['protomer_MW'] * echo_df['n_chains']
        echo_df['MW'] = echo_df['oligomer_MW']

    else:
        echo_df['MW'] = echo_df['protomer_MW'] 

    return echo_df

def add_vector_and_instrument_to_df(
    dataframe, 
    path_to_vector_folder, 
    sec_instrument,
    oligomers,
    culture_vol,
    culture_time,
    culture_temp,
    sec_column,
    sec_buffer,
    sec_buffer_pH,
):
    # Clone.
    dataframe["plasmid_seq"] = dataframe.apply(
        clone,
        path_to_vectors=path_to_vector_folder,
        axis=1,
    )
    dataframe["ORF"] = dataframe.apply(coding_seq, axis=1)

    # Protparam.
    dataframe = protparam(dataframe, oligomers, sec_buffer_pH)

    # Expression parameters.
    dataframe["culture_vol"] = float(culture_vol.split("x")[0]) * float(
        culture_vol.split("x")[1]
    )
    dataframe["culture_time"] = culture_time
    dataframe["culture_temp"] = culture_temp

    # Purification parameters.
    dataframe["sec_instrument"] = sec_instrument
    dataframe["sec_column"] = sec_column
    dataframe["sec_buffer"] = sec_buffer

    # OPTIONAL: define a category for later comparisons (e.g. designs made using different methods).
    dataframe["category"] = dataframe["Destination Plate Name"].astype(str)

#################
# SEC stuff
#################

def parse_chromatograms(
    sec_traces, 
    instrument=None,
):
    '''
    Parse SEC traces and return dataframe.
    
    sec_traces: list of paths to exported chromatograms from either AKTA or HPLC.
    
    instrument: the instrument used for SEC. Either 'akta', or 'hplc'.
    '''

    if instrument == 'akta':

        if os.path.isdir(sec_traces[0]):
            print(f'It looks like you used the HPLC, but selected sec_instrument = \'{instrument}\'.')
        
        # Extract wavelengths used for monitoring elutions
        fdata = pc_uni6(sec_traces[0]) # assumes the same wavelengths were used for all
        fdata.load()
        fdata.xml_parse()
        fdata.clean_up()
        Abs = [k for k in fdata.keys() if k[:2]=='UV']

        vol_data = {'vol' + wl.split('_')[-1]:[] for wl in Abs}
        Abs_data = {'A' + wl.split('_')[-1]:[] for wl in Abs}
        other_data = {'trace_id':[], 'fractions':[], 'frac_vol':[], 'frac_start_stop':[], 'flow_rate':[]}
        sec_data = {**vol_data, **Abs_data, **other_data}

        for trace in sec_traces:
            
            # Parse chromatograms.
            fdata = pc_uni6(trace)
            fdata.load()
            fdata.xml_parse()
            fdata.clean_up()
            
            # Fraction data.
            v_frac, frac = np.array(fdata['Fractions']['data']).T
            frac_id = [f for f in frac if '.' in f]
            frac_start_idx = np.where(np.isin(frac, frac_id))[0]
            frac_start = v_frac[frac_start_idx].astype(float)
            frac_stop = v_frac[frac_start_idx + 1].astype(float)
            frac_vol = np.subtract(frac_stop, frac_start)
    
            sec_data['trace_id'].append(int(trace.split()[-1].replace('.zip', '')))
            sec_data['fractions'].append(frac_id)
            sec_data['frac_vol'].append(frac_vol)
            sec_data['frac_start_stop'].append(list(zip(frac_start, frac_stop)))   
            sec_data['flow_rate'].append(np.array(fdata['System flow']['data']).T[1].mean())
            
            for A in Abs:
                wl = A.split('_')[-1]
                vol, absorbance = np.array(fdata[A]['data']).T # volume vs absorbance
                sec_data['A' + wl].append(absorbance)
                sec_data['vol' + wl].append(vol)


    elif instrument == 'hplc':
        
        # Parse fraction file.
        
        frac_file = [f for f in sec_traces if '.csv' in f]
        
        if len(frac_file) == 0:
            print(f'Fraction file not found (\'{instrument}\' selected). This file should be placed in the folder that contains the trace-specific sub-folders.')
        
        fractions = {}
        Abs = {} # save all wavelengths
        record = False
        with open(frac_file[0], 'r') as rf:

            for line in rf.readlines():
                line = line.strip()
                line = line.strip("\'\"")

                if line.startswith('Sample Vial Number,'):
                    # Sample Vial Number,D1F-A1,Flow,0.650 mL/min,Injection Volume,50.000,Injection Volume Unit,µL
                    _, trace_id, _ , flow_rate, _, injection_vol, _, _ = line.split(',')

                if line.startswith('Yes,Signal'):
                    # Yes,Signal A,280 nm,1 nm,No
                    _, DAD, _, _ = line.split()
                    signal, wavelength = DAD.split(',')
                    Abs[signal] = wavelength

                if line.startswith("#,Location") :
                    # ,Location,Start Time (min),End Time (min),Volume (mL),Trigger Reason
                    record = True
                    fraclist = []

                if line.startswith('Reference Detector ='):
                    # Reference Detector = DAD1 (Delay Time = 0.2180)
                    record = False
                    fraclist.append(float(flow_rate.split()[0]))
                    fraclist.append(float(injection_vol))
                    fractions[trace_id] = fraclist # dump fraction info

                if record:
                    fraclist.append(line.split(','))


        # Parse chromatograms.
        vol_data = {'vol'+v:[] for v in Abs.values()}
        Abs_data = {'A'+v:[] for v in Abs.values()}
        other_data = {'trace_id':[], 'fractions':[], 'frac_vol':[], 'frac_start_stop':[], 'flow_rate':[]}
        sec_data = {**vol_data, **Abs_data, **other_data}

        for folder in [f for f in sec_traces if '.csv' not in f]:
            trace_id = folder.split('_')[-1].split("/")[-1]
            trace_info = fractions[trace_id]
            flow_rate = trace_info[-2]

            # Fractions data.
            frac_info = np.array(trace_info[1:-2])
            frac_id = frac_info[:, 1]
            frac_start = frac_info[:, 2].astype(float) * flow_rate
            frac_stop = frac_info[:, 3].astype(float) * flow_rate
            frac_vol = np.subtract(frac_stop, frac_start)

            sec_data['trace_id'].append(trace_id)
            sec_data['fractions'].append(frac_id)
            sec_data['frac_vol'].append(frac_vol)
            sec_data['frac_start_stop'].append(list(zip(frac_start, frac_stop)))   
            sec_data['flow_rate'].append(flow_rate)

            for trace in glob.glob(folder + '/*.CSV'): # extract each wavelength
                #this is hacky and hardcoded. figure out a way to get A, B, etc w/o assumptions
                #print(trace)
                if 'DAD1' in trace:
                    wl = Abs[trace.split()[-1].split('DAD1')[1].replace('.CSV', '')]
                elif 'MWD1' in trace:
                    wl = Abs[trace.split()[-1].split('MWD1')[1].replace('.CSV', '')]
                elif 'FLD1' in trace:
                    wl = Abs[trace.split()[-1].split('FLD1')[1].replace('.CSV', '')]
                # Chromatogram data.
                chroma = pd.read_csv(trace, names=['time', 'Abs'])
                chroma['vol'] = chroma['time'] * flow_rate # convert time-base to volume-base

                sec_data['A' + wl].append(chroma['Abs'].to_numpy())
                sec_data['vol' + wl].append(chroma['vol'].to_numpy())
                
    else:
        print('Instrument name not recognised. Must be one of {akta, hplc}.')
        

    return pd.DataFrame.from_dict(sec_data)



def process_sec_data(
    sec_df, 
    wavelength=280,
):
    '''
    Process SEC chromatograms for a given wavelength.

    sec_df: dataframe of extracted SEC data (arrays).
    
    wavelength: wavelength used for the analysis (e.g. peak finding, etc..)
    '''

    processed_data = {
        'trace_id':[],
        'vol':[],
        'Abs':[],
        'vol_clipped':[],
        'Abs_corr':[],
        'Abs_norm':[],
        'main_peak':[],
        'main_peak_height':[],
        'tot_integral':[],
        'nPeaks':[],
        'monodisperse':[],
        'peaks':[],
        'peak_heights':[],
    }
    
    for idx, r in sec_df.set_index('trace_id').iterrows():
        vol = r['vol' + str(wavelength)]
        Abs = r['A' + str(wavelength)]
        
        processed_data['vol'].append(vol)
        processed_data['Abs'].append(Abs)
        processed_data['trace_id'].append(idx)

        # Integral and peaks.
        vmin, vmax = (0.9, 2.7) if np.max(vol) < 5 else (7, 21)
        idx_range = np.logical_and(vol>vmin, vol<vmax)
        vol_r = vol[idx_range]
        Abs_r = Abs[idx_range]
        Abs_corr = Abs_r - Abs_r.min()
        Abs_norm = Abs_corr / (Abs_r.max() - Abs_r.min())
        
        peak_idx = []
        peak_promi = 0.5
        height = 10
        for i in range(10):
            if len(peak_idx) == 0: # reduce prominence parameter in case of 'flat' traces.
                peak_idx, peak_props = signal.find_peaks(Abs_corr, height, prominence=peak_promi)
                peak_promi /= 2
            else:
                break
        for i in range(10):
            if len(peak_idx) == 0: # reduce prominence parameter in case of 'flat' traces.
                peak_idx, peak_props = signal.find_peaks(Abs_corr, height, prominence=peak_promi)
                height -= 1
            else:
                break
            
        main_peak = vol_r[peak_idx[np.argmax(peak_props['peak_heights'])]]
        main_peak_height = Abs_r[peak_idx[np.argmax(peak_props['peak_heights'])]]
        tot_integral = np.trapz(Abs_corr, x=vol_r)

        processed_data['vol_clipped'].append(vol_r)
        processed_data['Abs_corr'].append(Abs_corr)
        processed_data['Abs_norm'].append(Abs_norm)
        processed_data['main_peak'].append(main_peak)
        processed_data['main_peak_height'].append(main_peak_height)
        processed_data['tot_integral'].append(tot_integral)
        processed_data['nPeaks'].append(len(peak_idx))
        processed_data['monodisperse'].append(True if len(peak_idx)==1 else False)
        processed_data['peaks'].append(vol_r[peak_idx])
        processed_data['peak_heights'].append(peak_props['peak_heights'])

    processed_df = pd.DataFrame.from_dict(processed_data)

    # Reduce the number of trace points to allow inter-trace comparisons.
    vol = np.linspace(vmin, vmax, 100) 
    processed_df['light_idx'] = processed_df['vol_clipped'].apply(lambda x: [np.argmin(np.abs(x-v)) for v in vol])

    def sel_vol(r):
        return list(r['vol_clipped'][r['light_idx']])

    def sel_Abs(r):
        return list(r['Abs_corr'][r['light_idx']])

    processed_df['vol_light'] = processed_df.apply(sel_vol, axis=1)
    processed_df['Abs_light'] = processed_df.apply(sel_Abs, axis=1)
    processed_df['Abs_norm_light'] = processed_df['Abs_light']\
                    .apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))

    combined_df = sec_df.merge(processed_df, left_on='trace_id', right_on='trace_id')

    return combined_df

def link_traces_to_wells(echo_df, sec_data, sec_instrument):
    sec_df = pd.DataFrame()
    for run, sec_d in sec_data.items():

        w96 = [
            str(p) + "_" + r + str(c)
            for p in range(1, 10)
            for r in "ABCDEFGH"
            for c in range(1, 13)
        ]

        if sec_instrument == "akta":

            akta2well = {i + 1: w for i, w in enumerate(w96)}
            sec_d["well_id"] = sec_d["trace_id"].map(akta2well)
            sec_d["Destination Plate Name"] = sec_d["well_id"].apply(
                lambda x: int(x.split("_")[0])
            )
            sec_d["Destination Well"] = sec_d["well_id"].apply(lambda x: x.split("_")[1])

        if sec_instrument == "hplc":

            drawer2plate = {
                "D1F": 1,
                "D1B": 2,
                "D2F": 3,
                "D2B": 4,
                "D3F": 5,
                "D3B": 6,
                "D4F": 7,
                "D4B": 8,
            }
            sec_d["Destination Plate Name"] = sec_d["trace_id"].apply(
                lambda x: drawer2plate[x.split("-")[0]]
            )
            sec_d["Destination Well"] = sec_d["trace_id"].apply(lambda x: x.split("-")[1])
            sec_d["well_id"] = (
                sec_d["Destination Plate Name"].astype(str)
                + "_"
                + sec_d["Destination Well"]
            )

        sec_df = pd.concat([sec_df, sec_d])

    df = echo_df.merge(
        sec_df,
        left_on=["Destination Plate Name", "Destination Well"],
        right_on=["Destination Plate Name", "Destination Well"],
    ).reset_index(drop=True)
    return df

############################
# SEC calibration stuff
############################

def CI(
    x_values,
    y_values,
    xs,
    interval=0.95
):
    '''
    Calculate confidence interval of a linear fit to y_values vs. x_values.
    
    x_values: the x values to be fit.
    
    y_values: the corresponding y_values.
    
    xs: the x-range of the CI returned by the function.
    
    interval: the confidence interval (higher means less stringent).
    
    Returns: y_CI_lower, y_CI_upper
    
    CI calculation addapted from: from https://newbedev.com/show-confidence-limits-and-prediction-limits-in-scatter-plot
    '''
        
    x = x_values
    y = y_values

    slope, intercept = np.polyfit(x_values, y_values, 1)  # linear model adjustment

    y_model = np.polyval([slope, intercept], x_values)   # modeling...

    x_mean = np.mean(x_values)
    y_mean = np.mean(y_values)
    n = len(x_values)                 # number of samples
    m = 2                             # number of parameters
    dof = n - m                       # degrees of freedom
    t = stats.t.ppf((1+interval)/2, dof)       # Students statistic of interval confidence
    # (the confidence interval represents the probability that a parameter value falls within it.
    # A higher confidence interval is therefore less stringent, i.e. you're casting the net wider)

    residual = y_values - y_model

    std_error = (np.sum(residual**2) / dof)**.5   # standard deviation of the error

    # Confidence interval.
    ci = t * std_error * (1/n + (xs - x_mean)**2 / np.sum((x_values - x_mean)**2))**.5
    y_line = np.polyval([slope, intercept], xs)
    y_CI_upper = y_line + ci
    y_CI_lower = y_line - ci
    
    return y_CI_lower, y_CI_upper


def withinCI(
    r, 
    x, 
    y_lower, 
    y_upper
):
    '''Checks if main peak is within CI of the calibration curve'''

    log10mw = np.log10(r['MW'])
    Kav = r['main_peak_norm_retention']

    closest = np.argmin(np.abs(x - log10mw))
    up = y_upper[closest]
    low = y_lower[closest]

    if Kav <= up and Kav >= low:
        return True

    else:
        return False

def Vel2MW(Vel, sec_cal):
    Vo = sec_cal['Vo']
    Vc = sec_cal['Vc']
    intercept = sec_cal['intercept']
    slope = sec_cal['slope'] 
    return 10**((((Vel - Vo) / (Vc - Vo)) - intercept) / slope)
    
def calibrated_results(
    df, 
    sec_cal
):
    '''
    Compute calibrated result parameters using the information from a column calibration.
    
    df: the dataframe containing uncalibrated results. Calibrated result parameters will be appended to it. 
    
    sec_cal: generated from .json file containing the column calibration information.
    '''
    # Get SEC column calibration parameters.
    Vo = sec_cal['Vo']
    Vc = sec_cal['Vc']
    intercept = sec_cal['intercept']
    slope = sec_cal['slope']   
    all_mws = np.hstack([sec_cal['log10mw'], np.log10(df['MW'].to_numpy())])
    xs = np.linspace(all_mws.min()-0.2, all_mws.max()+0.2, 100)
    # cal_standards = Vo, Vc, intercept, slope, all_mws, xs

    # Compute parameters based on SEC calibration data
    df['expected_norm_retention'] = df['MW'].apply(
        lambda x: intercept + slope * np.log10(x)
    )

    df['expected_Vel'] = df['MW'].apply(
        lambda x: (Vc - Vo) * (intercept + slope * np.log10(x)) + Vo
    )

    df['main_peak_norm_retention'] = df['main_peak'].apply(
        lambda x: (x - Vo) / (Vc - Vo)
    )

    df['main_peak_estimated_MW'] = df['main_peak'].apply(
        lambda x: 10**((((x - Vo) / (Vc - Vo)) - intercept) / slope)
    )

    df['main_peak_agg_state'] = round(df['main_peak_estimated_MW'] / df['protomer_MW']).astype(int)

    # 95% confidence interval for correct retention volume.
    y_CI95_lower, y_CI95_upper = CI(sec_cal['log10mw'], sec_cal['Kav'], xs,interval=0.95)
    df['correct_Vel_95CI'] = df.apply(withinCI, x=xs, y_lower=y_CI95_lower, y_upper=y_CI95_upper, axis=1)

    # 99% confidence interval for correct retention volume.
    y_CI99_lower, y_CI99_upper = CI(sec_cal['log10mw'], sec_cal['Kav'], xs, interval=0.99)
    df['correct_Vel_99CI'] = df.apply(withinCI, x=xs, y_lower=y_CI99_lower, y_upper=y_CI99_upper, axis=1)
    
    return df, (y_CI95_lower, y_CI95_upper), (y_CI99_lower, y_CI99_upper)
  
    
def yields(df):
    '''
    Calculate protein yields from SEC chromatograms. 
    
    df: input dataframe to which the data gets appended.
    '''
    
    # Assign opitcal pathlength based on which instrument was used for SEC.
    df['path_length'] = df['sec_instrument'].apply(lambda x: 1 if x=='hplc' else 0.2)

    # Calculate yields.
    df['tot_yield'] = (df['tot_integral'] * 1e-3) / (df['path_length'] * df['OD280'])
    df['yield_per_Leq'] = (df['tot_yield'] * 1000) / df['culture_vol']

    if 'A205' in df.columns:
        df['tot_yield_from_205nm'] = (df['tot_integral'] * 1e-3) / (df['path_length'] * df['OD205'])
        df['yield_per_Leq_from_205nm'] = (df['tot_yield'] * 1000) / df['culture_vol'] 
        
    return df

def plot_all_traces(df):
    for cat, sel in df.groupby("category"):

        sel.reset_index(inplace=True)

        # Cluster chromatograms.
        vol_light = np.vstack(sel.loc[:, "vol_light"].to_numpy())
        Abs_light = np.vstack(sel.loc[:, "Abs_light"].to_numpy())
        clustered_idx = hierarchy.dendrogram(
            hierarchy.linkage(Abs_light, method="average", optimal_ordering=True),
            no_plot=True,
        )[
            "leaves"
        ]  # 'leaves' are the re-ordered indices

        # Plot.
        delta = np.max([np.max(x) for x in sel.Abs.to_numpy()]) / 20
        fig, ax = plt.subplots(ncols=2, figsize=(10, 5))
        for i, r in sel.iterrows():
            ax[0].plot(r.vol, r.Abs, color="C0", alpha=0.1)
            ax[1].fill_between(
                x=vol_light[clustered_idx[i]],
                y1=Abs_light[clustered_idx[i]] + i * delta,
                y2=i * delta,
                color="C0",
                alpha=0.1,
                zorder=i,
            )

        ax[0].set(xlabel="Retention vol. / mL",  ylabel=f"A280 / mAU")
        ax[1].set(xlabel="Retention vol. / mL", yticks=[])
        ax[1].spines["left"].set_visible(False)
        ax[0].set(title=f"{cat} ($N = {len(sel)}$)")
        plt.tight_layout()
        #plt.savefig(f"all_sec_traces_{cat}.png", dpi=300)
        plt.show()
        
def plot_calibrated_results(df, sec_calibration, CI95):
    with open(sec_calibration, "r") as f:
        sec_cal = json.load(f)

    Vo, Vc = sec_cal["Vo"], sec_cal["Vc"]
    intercept, slope = sec_cal["intercept"], sec_cal["slope"]
    all_mws = np.hstack([sec_cal["log10mw"], np.log10(df["MW"].to_numpy())])
    xs = np.linspace(all_mws.min() - 0.2, all_mws.max() + 0.2, 100)
    
    for c, sel in df.groupby("category"):

        # Plot calibration curve and CI.
        plt.scatter(
            [10**x for x in sec_cal["log10mw"]],
            sec_cal["Kav"],
            edgecolor="grey",
            linewidth=1,
            color="None",
            alpha=1,
        )

        plt.plot(
            10**xs,
            intercept + slope * xs,
            color="grey",
            linewidth=1,
        )

        plt.fill_between(
            10**xs,
            CI95[0],
            CI95[1],
            color="grey",
            alpha=0.2,
            edgecolor=None,
        )

        # Plot data.
        norm = plt.Normalize(sel.tot_yield.min(), sel.tot_yield.max())
        sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)

        ax = sns.scatterplot(
            data=sel,
            x="MW",
            y="main_peak_norm_retention",
            hue="tot_yield",
            palette="Spectral_r",
            zorder=10,
            edgecolor="k",
            linewidth=1,
            alpha=0.7,
        )

        plt.xscale("log")
        plt.xlabel("MW / Da")
        plt.ylabel("Norm. retention / a.u.")

        ax.get_legend().remove()
        plt.colorbar(sm, label="Tot. yield / mg",ax=ax)
        plt.title(c)
        plt.tight_layout()
        #plt.savefig(f"{c}_against_calibration.png", dpi=300)
        plt.show()

def plot_total_soluble_yield(df):
    fig, ax = plt.subplots(ncols=2, figsize=(10, 4))

    sns.histplot(
        data=df,
        x="tot_yield",
        stat="count",
        ax=ax[0],
        log_scale=True,
        element="step",
        alpha=0.1,
    )

    ax0 = ax[0].twinx()
    ax0.vlines(df.tot_yield.median(), 0, 0.5, color="salmon")
    ax0.scatter(df.tot_yield.median(), 0.5, color="salmon", zorder=10)
    sns.ecdfplot(data=df, x="tot_yield", alpha=0.7, ax=ax0)
    ax[0].spines["right"].set_visible(True)
    ax[0].set_xlabel("Total soluble yield / mg")
    ax[0].set_title(f'Median = {df["tot_yield"].median():.2f} mg\n($N$ = {len(df)})')

    sns.histplot(
        data=df,
        x="tot_yield",
        hue="category",
        stat="count",
        ax=ax[1],
        log_scale=True,
        element="step",
        alpha=0.1,
        legend=False,
    )

    sns.ecdfplot(
        data=df,
        x="tot_yield",
        hue="category",
        ax=ax[1].twinx(),
        legend=True,
        alpha=0.7,
    )

    ax[1].spines["right"].set_visible(True)
    ax[1].set_xlabel("Total soluble yield / mg")

    plt.tight_layout()
    #plt.savefig("tot_soluble_yields.png", dpi=300)
    plt.show() 
    
def plot_peak_distribution(df):
    xs = np.hstack(df.peaks.to_list())
    ys = np.hstack(df.peak_heights.to_list())
    g = sns.jointplot(
        x=xs,
        y=ys,
        kind="kde",
        ratio=4,
    )
    g.ax_joint.scatter(xs, ys, facecolors="none", edgecolor="salmon")
    g.ax_joint.set_xlabel("Peak retention volume / mL")
    g.ax_joint.set_ylabel("Normalized peak height")
    plt.tight_layout()
    #plt.savefig("peaks_distribution.png", dpi=300)
    plt.show()

def plot_peak_retention_over_norm_retention(df):
    norm = plt.Normalize(df.tot_yield.min(), df.tot_yield.max())
    sm = plt.cm.ScalarMappable(cmap="Spectral_r", norm=norm)

    sns.boxplot(
        data=df,
        x="category",
        y="expected_norm_retention",
        zorder=0,
        medianprops={"visible": True, "color": "lightgrey", "ls": "-", "lw": 10},
        whiskerprops={"visible": False},
        showfliers=False,
        showbox=False,
        showcaps=False,
    )

    ax = sns.swarmplot(
        data=df,
        x="category",
        y="main_peak_norm_retention",
        hue="tot_yield",
        palette="Spectral_r",
        alpha=0.7,
        edgecolor="k",
        linewidth=1,
        size=6,
    )

    ax.get_legend().remove()
    ax.figure.colorbar(sm, label="Tot. yield / mg")
    ax.set(ylim=[-0.1, 1], xlabel="", ylabel="Norm. retention / a.u.")
    plt.tight_layout()
    #plt.savefig("norm_retention_by_category.png", dpi=300)
    plt.show()

def plot_aggregation_state_by_plate(df):
    sns.stripplot(
        data=df.sort_values(by="main_peak_agg_state"),
        x="main_peak_agg_state",
        y="main_peak_height",
        hue="category",
        alpha=0.5,
        edgecolor="k",
        linewidth=1,
        size=6,
        dodge=True,
    )

    plt.xlabel("Aggregation state")
    plt.ylabel("Peak height / mAU")
    plt.tight_layout()
    #plt.savefig("aggregation_states.png", dpi=300)
    plt.show()
    
def plot_aggregation_state_by_well(df, sec_cal):
    for plate, sel in df.groupby("Destination Plate Name"):

        states = [str(i) for i in range(9)] + ["9+"]

        agg_df = pd.DataFrame(
            index=sel["Destination Well"].to_list(),
            columns=states,
        )

        agg_df = agg_df.merge(
            sel[["tot_yield", "Destination Well"]],
            left_index=True,
            right_on="Destination Well",
        ).set_index("Destination Well", drop=True)
        
        with pd.option_context("future.no_silent_downcasting", True):
            agg_df = agg_df.fillna(0).astype(np.float64)
        
        for i, r in sel.iterrows():
            agg_states = [round(Vel2MW(v, sec_cal) / r.protomer_MW) for v in r.peaks]
            proportions = r.tot_yield * r.peak_heights / r.peak_heights.sum()
            
            for agg_st, prop in zip(agg_states, proportions):

                if agg_st < 9:
                    agg_df.loc[r["Destination Well"], str(agg_st)] += prop.astype(np.float64)

                else:
                    agg_df.loc[r["Destination Well"], "9+"] += prop.astype(np.float64)
        agg_df.sort_values(
            by=["tot_yield"] + states,
            ascending=False,
            inplace=True,
        )

        n_subplots = np.ceil(len(agg_df) / 48).astype(int)
        fig, ax = plt.subplots(
            nrows=n_subplots,
            figsize=(len(agg_df) / n_subplots * 0.3, 5 * n_subplots),
            sharey=True,
        )

        if n_subplots == 1:
            ax = [ax]

        for i in range(n_subplots):
            agg_df.iloc[48 * i : 48 * (i + 1)][states].plot(
                ax=ax[i],
                kind="bar",
                stacked=True,
                edgecolor="k",
                alpha=0.80,
                cmap="magma",
            )

            ax[i].legend(bbox_to_anchor=(1, 1), title="Agg. state")
            ax[i].set(xlabel="", ylabel="Total soluble yield / mg")
            ax[i].set_title(f"Plate {plate} ({48*i+1}-{48*(i+1)})")
        plt.tight_layout()
        #plt.savefig(f"plate_{plate}_agg_per_well.png", dpi=300)
        plt.show()

########################
# SEC pooling stuff
########################
    
def select_fractions(
    w,
    r,
    pooled_df,
    manual_edits,
    sec_cal,
    wl,
    n_fractions=1,
    how='nearest',
    adjacency=0.025,
    wl2=None
):
    '''
    Automatically select fraction(s) for pooling.
    
    w: the well id for the chromatogram.
    
    r: the row of the dataframe corresponding to that well id.
    
    pooled_df: dataframe for saving the data.
    
    manual_edits: array of manual edits for fraction picking. Overrides automatic picking for these traces.
    
    sec_cal: generated from .jason containing calibration information.
    
    wl: the chromatogram wavelength.
    
    n_fractions: the number of adjacent fractions to pick.
    
    how: {'largest', 'nearest'} automatic fraction picking logic. 
          - 'largest' will pick the fraction(s) centered around the fraction with the largest integral. 
          - 'nearest' will pick the fractions(s) centered around the fraction nearest to the expected retention volume.
    '''
    # Extract relevant data from that trace.
    vol = r.vol_clipped
    Abs = r.Abs_corr
    frac = r.fractions
    vol_start_stop = r.frac_start_stop
    peaks_vol = r.peaks
 
    if len(frac) == 0:
        pooled_df.loc[w, 'pooled_fractions'] = []
        pooled_df.loc[w, 'pooled_frac_vol'] = []
        pooled_df.loc[w, 'pooled_vol'] = 0
        pooled_df.loc[w, 'pooled_integral'] = 0
        pooled_df.loc[w, 'pooled_peak_Abs'] = 0
        pooled_df.loc[w, 'pooled_peak_vol'] = 0

    else:
        # Compute properties for each fraction.
        fractions = []
        frac_edges = []
        frac_integrals = []
        frac_peaks_Abs = []
        frac_peaks_vol = []
        for i, (v1, v2) in enumerate(vol_start_stop):
            i_range = np.argwhere(np.logical_and(vol>=v1, vol<=v2)).flatten()

            if len(i_range) > 0:
                fractions.append(frac[i])
                frac_edges.append([v1, v2])
                frac_integrals.append(np.trapz(Abs[i_range], x=vol[i_range]))
                p, _ = signal.find_peaks(Abs[i_range], prominence=0.01)
                
                if len(p) > 0:
                    frac_peaks_Abs.append(Abs[i_range][p].max())
                    frac_peaks_vol.append(vol[i_range][p[np.argmax(Abs[i_range][p])]])
                
                else:
                    frac_peaks_Abs.append(0)
                    frac_peaks_vol.append(0)
                
            else:
                pass
        if fractions == []: #if no fractions are picked return nothing
            return None
                    
            
        if how=='largest': # pick fraction with the largest integral
            picks = np.isin(frac_integrals, np.max(frac_integrals))
        

        elif how=='nearest': # pick fraction containin the peak nearest to the expected retention volume
            nearest_peak = peaks_vol[np.abs(peaks_vol - r.expected_Vel).argmin()]
            deltas = np.abs(np.subtract(frac_edges, nearest_peak)).min(axis=1)
            nearest_frac_idx = np.where(deltas==deltas.min())[0]
            
            if len(nearest_frac_idx) > 1: # in case of equidistance
                picks = np.array([], dtype=bool)
                
                for v in frac_edges:
                    picks = np.concatenate([picks, [np.logical_and(nearest_peak>=v[0], nearest_peak<=v[1])]])
                
            else:

                picks = np.zeros(len(fractions), dtype=bool)
                picks[nearest_frac_idx] = True           

        else:
            print('Automatic picking option not recognized. Must be either \'largest\' or \'nearest\'.')
            
        
        # Pick adjacent fraction(s) with the largest integrals.
        delta_vol_tol = adjacency # volume tolerance (in mL) for counting the previous/next fraction adjacent.   
        for i in range(n_fractions-1):

            left_integral, right_integral = 0, 0

            if i >= (len(frac_integrals) - 1):
                pass

            else:
                
                leftmost_idx = np.where(picks)[0].min()
                
                if leftmost_idx > 0 :
                    if (frac_edges[leftmost_idx - 1][1] + delta_vol_tol) >= frac_edges[leftmost_idx][0]:
                        left_integral = frac_integrals[leftmost_idx - 1]
                        
                    else: left_integral = 0

                else: left_integral = 0

                rightmost_idx = np.where(picks)[0].max()
                
                if rightmost_idx < (len(frac_integrals) - 1):
                    if (frac_edges[rightmost_idx][1] + delta_vol_tol) >= frac_edges[rightmost_idx + 1][0]:
                        right_integral = frac_integrals[rightmost_idx + 1]
                        
                    else: right_integral = 0

                else: right_integral = 0


            if left_integral + right_integral == 0:
                pass

            else:
                picks += np.isin(frac_integrals, np.max([left_integral, right_integral]))


        # Manual corrections, if any.
        if w in manual_edits[:, 0]:
            pooled_fractions = manual_edits[np.argwhere(manual_edits[:, 0]==w)].flatten()[1]
            picks = np.isin(fractions, pooled_fractions)

            
        pooled_frac_edges = np.array(frac_edges)[picks]
        picked_peaks = np.array(frac_peaks_Abs)[picks]
        pooled_peak_Abs = picked_peaks.max()
        pooled_peak_vol = np.array(frac_peaks_vol)[picks][picked_peaks.argmax()]

        
        # Save data.
        pooled_df.loc[w, 'pooled_fractions'] = np.array(fractions)[picks]
        pooled_df.loc[w, 'pooled_frac_vol'] = np.diff(pooled_frac_edges, axis=1).flatten()
        pooled_df.loc[w, 'pooled_vol'] = np.diff(pooled_frac_edges, axis=1).sum()
        pooled_df.loc[w, 'pooled_integral'] = np.array(frac_integrals)[picks].sum() 
        pooled_df.loc[w, 'pooled_peak_Abs'] = pooled_peak_Abs
        pooled_df.loc[w, 'pooled_peak_vol'] = pooled_peak_vol
        

        # Shade pooled fraction(s).
        for v1, v2 in pooled_frac_edges:
            i_range = np.argwhere(np.logical_and(vol>v1, vol<v2)).flatten()
            plt.fill_between(vol[i_range], Abs[i_range], color='C0', edgecolor='None', alpha=0.3)


        # Indicate fractions.
        plt.vlines(vol_start_stop, -0.2 * Abs.max(), 0, lw=0.5, color='grey')
        for i, frac_label in enumerate(fractions):
            
            if (frac_edges[i][0] >= vol.min()) & (frac_edges[i][1] <= vol.max()):
                plt.text(frac_edges[i][0], -0.2 * Abs.max(), frac_label,
                         fontsize=12, fontfamily='monospace', color='grey', rotation=45)

    
    # Plot.
    plt.plot(vol, Abs)
    if wl2 != None:
        vol2=r['vol'+str(wl2)]
        Abs2=r['A'+str(wl2)]
        plt.plot(vol2, Abs2) 
    
    # Indicate peaks.
    for peak in r.peaks:
        plt.text(peak, Abs[np.argwhere(vol==peak)], f'{Vel2MW(peak, sec_cal)*1e-3:.0f}',
                fontsize=12, fontfamily='monospace', color='grey', rotation=45,
                va='bottom', ha='center')

        
    # Indicate anticipated retention volume.
    plt.vlines(r.expected_Vel, 0, Abs.max(), lw=10, color='salmon', alpha=0.3, zorder=0)
    
    
    # Polish.
    name = r.Name if len(r.Name) < 30 else r.Name.replace(r.Name[30:], "...")
    plt.title(f'{w} /// {r.MW * 1e-3:.0f} kDa\n{name}' + r'$\in$' + f'{r.Vector}')
    plt.xlabel('Retention vol. / mL')
    plt.ylabel(f'A{wl} / mAU')
    plt.xlim([vol.min(), vol.max()])
    plt.ylim([Abs.min() - 0.25 * Abs.max(), Abs.max() + 0.2 * Abs.max()])
    plt.show()
    
    
    if len(frac) == 0:
        print(f'[!] No fraction collected for {w}')
        
    else:
        # Display next best fraction.
        remaining_integrals = np.array(frac_integrals)[~picks]
        if len(remaining_integrals) > 0:
            next_best_fraction = np.array(fractions)[np.isin(frac_integrals, remaining_integrals.max())]
            print(f"Next best fraction:\n['{w}', {next_best_fraction}],")

    return

#################
# CONCENTRATION
#################
def calculate_conc_and_add_to_database(df, sec_cal, expression_name):
    df["conc_uM"] = (
        1e3 * df["pooled_integral"] / (df["pooled_vol"] * df["e280"] * df["path_length"])
    )
    df["conc_mgmL"] = (
        1e-3 * df["pooled_integral"] / (df["pooled_vol"] * df["OD280"] * df["path_length"])
    )
    df["pooled_mg"] = df["conc_mgmL"] * df["pooled_vol"]
    df["pooled_peak_estimated_MW"] = df["pooled_peak_vol"].apply(Vel2MW, args=(sec_cal,))
    df["pooled_peak_agg_state"] = round(
        df["pooled_peak_estimated_MW"] / df["protomer_MW"]
    ).astype(int)
    df["pooled_is_main"] = df["main_peak"] == df["pooled_peak_vol"]
    user = getpass.getuser()
    file = f"/projects/eblocks/sec/{user}_{today}_{len(df)}_{expression_name}_eblock_sec.csv"
    df.to_csv(file)
    return df
    
def plot_pooled_vol_and_conc(df):
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10, 6))
    ax = ax.reshape(-1)

    for i, (k, v) in enumerate(
        {
            "pooled_vol": "Pooled volume / mL",
            "pooled_mg": "Pooled amount / mg",
            "conc_uM": "Concentration / $\mu$M",
            "conc_mgmL": "Concentration / mg mL$^{-1}$",
        }.items()
    ):

        sns.histplot(
            df[k],
            ax=ax[i],
            element="step",
            bins=15,
            alpha=0.1,
            log_scale=False if "vol" in k else True,
        )
        sns.ecdfplot(df[k], ax=ax[i].twinx(), alpha=0.5, color="salmon")
        ax[i].spines["right"].set_visible(True)
        ax[i].set(xlabel=v)

    plt.tight_layout()
    #plt.savefig("pooled_results.png", dpi=300)
    plt.show()
    
def correct_pooled(r, desired_conc, concentration_type = "uM", max_vol=1.5, max_pooled = 10000, normalize = False):
    '''
    Correct pooled volumes for samples that are too high in concentration 
    (pooled volumes are reset to ensure; buffer + sample = max_vol and concentration will be correct)
    
    desired_conc: the concentration in uM that you want to set
    max_vol: the maximum volume of the combined fractions and buffer that is acceptable (in mL).
    '''
    if (r.pooled_vol * (r[f'conc_{concentration_type}'] - desired_conc) / desired_conc) + r.pooled_vol > max_vol:
        #if the pooled volume is too high to normalize to desired concentration in destination plate
        if (max_vol * desired_conc)/r[f'conc_{concentration_type}'] < 0.01:
            #if it would have to be diluted down to less than 10 uL to reach the desired concentration
            return 0.01 # just return 10 uL. Maybe should print a warning since this triggering feels like someone messed up
        elif normalize:
            return min(((max_vol * desired_conc)/r[f'conc_{concentration_type}'])*.8,max_pooled/1000)
        else:
            return min((max_vol * desired_conc)/r[f'conc_{concentration_type}'],max_pooled/1000)
    elif normalize:
        return min((r.pooled_vol*.8), max_vol, (max_pooled/1000))
    else:
        return min((r.pooled_vol), max_vol, (max_pooled/1000))

def correct_buffer(r, max_vol=1.5):
    '''
    Correct buffer volumes for samples that are either:
    - too low in concentration (negative buffer volumes are reset to zero).
    - too high in concentration (large buffer volumes are reset to ensure; buffer + sample = max_vol)
    
    max_vol: the maximum volume of the combined fractions and buffer that is acceptable (in mL). Above that buffer_vol is modified.
    '''
    
    if r.buffer_vol < 0:
        return 0
        
    elif np.isinf(r.conc_uM):
        #print(f'Well {r.position} has infinite concentration. This is likely due to no absorbance at 280 nM. Not normalizing')
        return 0
    
    elif r.buffer_vol + r.pooled_vol_corr > max_vol:
        return max_vol - r.pooled_vol_corr
    
    else:
        return r.buffer_vol
def calculate_dilutions(df, desired_conc, concentration_type = "uM", max_volume=1.5, max_pooled_vol=10000, verbose=True, normalize = False):
    df['pooled_vol_corr'] = df.apply(correct_pooled, args=(desired_conc, concentration_type, max_volume,max_pooled_vol,normalize), axis=1)
    df['buffer_vol'] = df.pooled_vol_corr * (df[f'conc_{concentration_type}'] - desired_conc) / desired_conc
    df['buffer_vol'] = df.apply(correct_buffer, args=(max_volume,), axis=1)
    df[f'normed_conc_{concentration_type}'] = df[f'conc_{concentration_type}'] * df.pooled_vol_corr / (df.pooled_vol_corr + df.buffer_vol)

    frac_correct = np.isclose(
        df[f'normed_conc_{concentration_type}'].to_list(), desired_conc, atol=0.05
    ).sum() / len(df)
    if verbose:
        print(f"{100 * frac_correct:.1f}% of wells will have the desired concentration.")
        
def find_opt_concentration(dataframe, concentration_type = "uM"):
    '''     
    Find the normalisation concentration that maximises the number of wells that will be normalised to said concentration.
    '''
    print("This functionality was changed in version 8. Make sure you've grabbed the most recent copy of the template from /software/lab/cowboy/")
    df = dataframe.copy()
    desired_conc = np.logspace(-2, 2, 100) # uM
    frac_correct = []
    for d_conc in desired_conc:
        calculate_dilutions(df, d_conc, concentration_type, 1.5, verbose=False)
        frac_correct.append(np.isclose(df[f'normed_conc_{concentration_type}'].to_list(), d_conc, atol=0.05).sum()/ len(df))
    max_corr = frac_correct[np.argmax(frac_correct)]
    opt_conc = desired_conc[np.argmax(frac_correct)]
    percent_90 = 0
    percent_80 = 0
    for i, num in enumerate(frac_correct):
        if num == max_corr:
            if opt_conc < desired_conc[i]:
                opt_conc = desired_conc[i]
        if num >= .8:
            if desired_conc[percent_80] < desired_conc[i]:
                percent_80 = i
            if num >=.9:
                if desired_conc[percent_90] < desired_conc[i]:
                    percent_90 = i

    print(f'Choosing {opt_conc} {concentration_type} will give {max_corr*100}% correct wells. This is the highest percentage')
    print(f'Choosing {desired_conc[percent_90]:.2f} {concentration_type} will give {frac_correct[percent_90]*100}% correct wells.' )
    print(f'Choosing {desired_conc[percent_80]:.2f} {concentration_type} will give {frac_correct[percent_80]*100}% correct wells.')

def plot_conc_and_vol_after_norm(df, concentration_type = "uM"):
    fig, ax = plt.subplots(ncols=2, figsize=(10, 3.5))
    ax = ax.reshape(-1)

    for i, (k, v) in enumerate(
        {
            f"normed_conc_{concentration_type}": f"Normalised concentration / {concentration_type}",
            "buffer_vol": "Volume buffer / mL",
        }.items()
    ):

        sns.histplot(df[k], ax=ax[i], element="step", alpha=0.1, bins=30)
        sns.ecdfplot(df[k], ax=ax[i].twinx(), alpha=0.5, color="salmon")
        ax[i].spines["right"].set_visible(True)
        ax[i].set(xlabel=v)

    plt.tight_layout()
    #plt.savefig("normalization_results.png", dpi=300)
    plt.show()
#################
# OT-2 stuff
#################

def gen_ot2_script(
    df,
    sec_instrument,
    plate_type,
    normalize=True,
    separate_non_normed=False,
    filename="",
    path="",
    date=today,
    template="/net/software/lab/cowboy/SEC_pool_and_norm_v3.py",
):
    """
    Generate OT-2 script for pooling and normalising SEC fractions.

    df: dataframe containing the experimental data.

    filename: filename for the generated script.

    date: string of the date that gets append to the output files.

    template: path to an OT-2 protocol template.
    """
    
    user = getpass.getuser()
    file_name = f"{path}{date}_{filename}{user}_sec_pool_and_norm"
    if separate_non_normed == False:
        print(
            f"{file_name}.py contains transfers from {len(np.hstack(df.pooled_fractions.to_list()))} fractions to {len(df)} destination wells"
        )

    ot2_transfers = ""
    non_normed_transfers = ""
    for i, r in df.iterrows():
        if normalize:
            volume = r.pooled_vol_corr
            if r.pooled_fractions.size > 1:
                print("Warning! Currently normalization only works properly if you only pull a single fraction per run")
                return
        else:
            volume = r.pooled_frac_vol
        pooled_frac_list = " ".join(
            [
                f"'{a}',"
                for a in str(r.pooled_fractions)
                .replace("['", "")
                .replace("']", "")
                .split("' '")
                if a != ""
            ]
        )

        pooled_frac_vol_list = [
            float(a) * 1e3
            for a in str(volume).replace("[", "").replace("]", "").split(" ")
            if a != ""
        ]
        buffer_vol = r.buffer_vol * 1e3
        writestr = f"[ [{pooled_frac_list}], {pooled_frac_vol_list}, {r['Destination Plate Name']}, '{r['Destination Well']}', {buffer_vol:3.3f}],\n"
        if normalize and separate_non_normed and buffer_vol == 0.000:
            non_normed_transfers += writestr
        else:
            ot2_transfers += writestr
    if normalize and separate_non_normed:
        print(
            f"{file_name}.py contains transfers to {len(ot2_transfers.splitlines())} destination wells"
        )
    print(ot2_transfers)
    ot2_template = open(template, "r").readlines()

    with open(file_name + ".py", "w") as f:
        for l in ot2_template:
            modified_line = l.replace("### PROTOCOL_NAME ###", file_name)
            modified_line = modified_line.replace("### USER_NAME ###", user)
            modified_line = modified_line.replace("### INSTRUMENT ###", sec_instrument)
            modified_line = modified_line.replace("### PLATE TYPE ###", plate_type)
            modified_line = modified_line.replace(
                "### NORMALIZE ###", "True" if normalize else "False"
            )
            modified_line = modified_line.replace("### TRANSFERS ###", ot2_transfers)
            f.write(modified_line)
    if separate_non_normed:
        print(
            f"{file_name}_non_normed.py contains transfers to {len(non_normed_transfers.splitlines())} destination wells"
        )
        print(non_normed_transfers)
        with open(file_name + "_non_normed.py", "w") as f:
            for l in ot2_template:
                modified_line = l.replace("### PROTOCOL_NAME ###", file_name)
                modified_line = modified_line.replace("### USER_NAME ###", user) 
                modified_line = modified_line.replace(
                    "### INSTRUMENT ###", sec_instrument
                )
                modified_line = modified_line.replace("### NORMALIZE ###", "True")
                modified_line = modified_line.replace(
                    "### TRANSFERS ###", non_normed_transfers
                )
                f.write(modified_line)

    return

