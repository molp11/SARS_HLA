from os import listdir
from os.path import join
import numpy as np
import pandas as pd
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

def aff2norm(a):
    return 1 - np.log(a)/np.log(50000)

def norm2aff(n):
    return np.exp((1 - n) * np.log(50000))

def mhc_filter(ba, el, el_lim, ba_lim):
    output = ba.where(el > el_lim, 0.).where(ba > ba_lim, 0.)
    return output.where (output == 0., 1.)

def reduce_df(df):
    df.index = df.index.set_levels([str(x) for x in df.index.levels[0]], level=0)
    residues = df.index.get_level_values(0)
    unique = residues[residues.duplicated() == False]
    idx = pd.IndexSlice
    output = df.loc[idx[unique, :, :], idx[:, :, :, :]].droplevel(2).sort_index()
    duplicated = residues[residues.duplicated() == True].unique()
    dupDf = df.loc[duplicated]
    output2 = pd.DataFrame()
    for c in dupDf.columns:
        Series = dupDf[c].loc[idx[:, :, c[2].split(' ')+['WT']]]
        dupIndex = Series[Series.index.get_level_values(0).duplicated() == True].index.get_level_values(0)
        uniqueIndex = []
        for x in Series.index.get_level_values(0).unique():
            if x not in dupIndex:
                uniqueIndex.append(x)
        dupSeries = Series.loc[idx[dupIndex.values, :, :]]
        uniqueSeries = Series.loc[idx[uniqueIndex, :, :]]
        tmp = pd.DataFrame(dupSeries).loc[idx[:, :, [x for x in dupSeries.index.get_level_values(2) if x in c[2]]]].append(pd.DataFrame(uniqueSeries)).droplevel(2).sort_index()
        if len(output2) == 0:
            output2 = tmp
        else:
            output2[c] = tmp[c].copy()
    return output.append(output2)

def get_coverage(reduced):
    length = max([int(x.split(',')[1].strip(')')) for x in reduced.index.levels[0]])
    coverage_df = pd.DataFrame()
    srange = [ix[0].replace('(', '').replace(')', '').replace(',', '').split(' ') for ix in reduced.index]
    for s in reduced.columns:
        coverage = np.zeros(length)
        if np.max(reduced[s]) > 0:
            for r, i in zip(reduced[s], srange):
                start = int(i[0])
                stop = int(i[1])
                for x in  range(start, stop-1):
                    if coverage[x] == 0. and r == 1.:
                        coverage[x] = r
        coverage_df[s] = coverage
    return coverage_df

protein_names = ['NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'NSP6', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'NSP11', 
                 'NSP12', 'NSP13', 'NSP14', 'NSP15', 'NSP16', 'Spike', 'NS3', 'M', 'EE', 'NS6', 'NS7a', 'NS7b',
                 'NS8', 'NS0b', 'NS99c', 'N']

for ba_lim in np.linspace(0.5, 0.95, 10):
    for el_lim in np.linspace(0.5, 0.95, 10):
        
        el_lim = round(el_lim, 2)
        ba_lim = round(ba_lim, 2)
        
        # MHCI
        master = pd.DataFrame()
        for tq, protein in zip(tqdm(range(len(protein_names)), desc='Filtering  MHCI (EL: ' + str(el_lim) + ', BA: ' + str(ba_lim) + ')'), protein_names):
            ba_files = [x for x in listdir('output/netMHCpan/pkl') if x.split('_')[0] == protein and 'BA.pkl' in x]
            el_files =  [x.replace('BA.pkl', 'EL.pkl') for x in ba_files]
            proteinDf = pd.DataFrame()
            for fba, fel in zip(ba_files, el_files):
                ba = pd.read_pickle(join('output/netMHCpan/pkl', fba))
                el = pd.read_pickle(join('output/netMHCpan/pkl', fel))
                filtered = mhc_filter(ba, el, el_lim, ba_lim)
                if len(filtered) > 0:
                    reduced = reduce_df(filtered)
                    if np.max(reduced.values) == 0:
                        continue
                    coverage_df = get_coverage(reduced)
                    coverage_df.index = pd.MultiIndex.from_tuples([(protein, x) for x  in coverage_df.index])
                    coverage_df.columns = pd.MultiIndex.from_tuples(coverage_df.columns)
                    if len(proteinDf) == 0:
                        proteinDf = coverage_df
                    else:
                        proteinDf = proteinDf + coverage_df
            proteinDf = proteinDf.where(proteinDf == 0, 1)
            if len(master) == 0:
                master = proteinDf
            else:
                master = master.append(proteinDf)
        master.to_pickle(join('output/netMHCpan', 'MHCI_BA' + str(ba_lim*100) + '_EL' + str(ba_lim*100) + '.pkl'))

        # MHCII
        master = pd.DataFrame()
        for tq, protein in zip(tqdm(range(len(protein_names)), desc='Filtering MHCII (EL: ' + str(el_lim) + ', BA: ' + str(ba_lim) + ')'), protein_names):
            ba_files = [x for x in listdir('output/netMHCIIpan/pkl') if x.split('_')[0] == protein and 'BA.pkl' in x]
            el_files =  [x.replace('BA.pkl', 'EL.pkl') for x in ba_files]
            proteinDf = pd.DataFrame()
            for fba, fel in zip(ba_files, el_files):
                ba = pd.read_pickle(join('output/netMHCIIpan/pkl', fba))
                el = pd.read_pickle(join('output/netMHCIIpan/pkl', fel))
                filtered = mhc_filter(ba, el, el_lim, ba_lim)
                if len(filtered) > 0:
                    reduced = reduce_df(filtered)
                    if np.max(reduced.values) == 0:
                        continue
                    coverage_df = get_coverage(reduced)
                    coverage_df.index = pd.MultiIndex.from_tuples([(protein, x) for x  in coverage_df.index])
                    coverage_df.columns = pd.MultiIndex.from_tuples(coverage_df.columns)
                    if len(proteinDf) == 0:
                        proteinDf = coverage_df
                    else:
                        proteinDf = proteinDf + coverage_df
            proteinDf = proteinDf.where(proteinDf == 0, 1)
            if len(master) == 0:
                master = proteinDf
            else:
                master = master.append(proteinDf)
        master.to_pickle(join('output/netMHCIIpan', 'MHCII_BA' + str(ba_lim*100) + '_EL' + str(ba_lim*100) + '.pkl'))  