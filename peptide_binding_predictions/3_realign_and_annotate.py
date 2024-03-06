from os import listdir
from os.path import join, exists
import pickle
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

class nestedDict(dict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

def add_column(df, score_column, data, records, gene, protein, length, allele, mutations, strain, sample):
    '''
    Realign peptides to SARS-CoV-2 genome and add to multiindexed dataframe
    '''
    if len(data) > 0:
        tuples, scores = [], []
        for r in records:
            if strain in r.id.split('|')[3]:
                for d in data.iterrows():
                    seq = d[1].Peptide
                    sc = d[1][score_column]
                    i = -1
                    mut_names = ''
                    if seq in r.seq:
                        i = r.seq.find(seq)
                        mut = mutation_dict[r.id.split('|')[3]][protein]
                        for m in mut:
                            if m != 'WT':
                                if int(m[1:-1]) in range(i+1, i+length+1):
                                    mut_names = [m]
                    else:
                        for rx in records:
                            if seq in rx.seq and 'Wuhan' not in rx.id:
                                i = rx.seq.find(seq)
                                mut = mutation_dict[rx.id.split('|')[3]][protein]
                                for m in mut:
                                    if m != 'WT':
                                        if int(m[1:-1]) in range(i+1, i+length+1):
                                            mut_names = [m]
                    if len(mut_names) == 0:
                        mut_names = ['WT']
                    scores.append(sc)
                    tuples.append(((i, i+length), str(seq), str(' '.join(mut_names))))
        index = pd.MultiIndex.from_tuples(tuples, names = ['residues', 'peptide', 'mutation'])
        columns = pd.MultiIndex.from_tuples([(sample, strain, str(' '.join(mutations)), allele)], 
                                            names=['sample', 'strain', 'mutation', 'allele'])
        tmp = pd.DataFrame({columns[0] : scores}, index=index)
        df[columns[0]] = tmp[columns[0]]
        assert True not in df.isna().any().values
    
    return df
        

df = pd.read_csv(join('input', 'data.csv'))

fastas = [f for f in listdir(join('input', 'fasta'))] 
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
mutation_dict = nestedDict()
for tq, fasta in zip(tqdm(range(len(fastas)), desc='Identifying mutations:'), fastas):
    protein = fasta.split('.')[0]
    records = list(SeqIO.parse(join('input/fasta', fasta), 'fasta'))
    
    # find Wuhan strain in fasta
    wuhan = ''
    for r in records:
        if 'Wuhan' in r.id:
            wuhan = r
            
    # compare patient strains to original wuhan
    for r in records:
        if 'Wuhan' not in r.id:
            gisaid_epi_isl = r.id.split('|')[3]
            strain = '/'.join(r.id.split('|')[1].split('/')[1:])
            for ix in range(min([len(r.seq), len(wuhan.seq)])):
                s1, s2 = wuhan.seq, r.seq
                if s1[ix] in AA and s2[ix] in AA and s1[ix] != s2[ix]:
                    mut = s1[ix] + str(ix + 1) + s2[ix]
                    if protein not in mutation_dict[gisaid_epi_isl]: 
                        mutation_dict[gisaid_epi_isl][protein] = [mut]
                    elif mut not in mutation_dict[gisaid_epi_isl][protein]:
                        mutation_dict[gisaid_epi_isl][protein].append(mut)   
                    if protein not in mutation_dict[strain]: 
                        mutation_dict[strain][protein] = [mut]                    
                    elif mut not in mutation_dict[gisaid_epi_isl][protein]:
                        mutation_dict[strain][protein].append(mut)
            if len(mutation_dict[gisaid_epi_isl][protein]) == 0:
                mutation_dict[gisaid_epi_isl][protein] = ['WT']
            if len(mutation_dict[strain][protein]):
                mutation_dict[strain][protein] = ['WT']
                
with open(join('input', 'UM_mutation_dict.pkl'), 'wb') as outfile:
    pickle.dump(mutation_dict, outfile)
    
print('Realigning and annotating netMHCpan predictions...')
NetMHCpan = pd.read_csv(join('output/netMHCpan/', 'NetMHCpan-output.csv'))
for length in range(8, 13):
    for tq, fasta in  zip(tqdm(range(len(fastas)), desc='MHCI ' + str(length) + 'mers'),  fastas):
        protein = fasta.split('.')[0]
        proteinData = NetMHCpan[NetMHCpan.Protein == protein]
        records = list(SeqIO.parse(join('input/fasta', fasta), 'fasta'))
        df_BA, df_EL =  pd.DataFrame(), pd.DataFrame()
        for row in df.iterrows():
            strain = row[1].gisaid_epi_isl
            sample = row[1]['SAMPLE_ID']
            mutations = mutation_dict[strain][protein]
            for col in ['A-1', 'A-2', 'B-1', 'B-2', 'C-1', 'C-2']:
                gene = col.split('-')[0]
                allele = row[1][col].strip(' ').split('(')[0]
                allele_converted = 'HLA-' + ':'.join(allele.replace('*', '').split(':')[:2])
                data = proteinData[(proteinData.Peptide_length == length) & (proteinData.HLA_type.str.contains(allele_converted))]
                if len(data) == 0 and allele_converted in NetMHCpan.HLA_type.unique():
                    data = NetMHCpan[(NetMHCpan.Peptide_length == length) & (NetMHCpan.HLA_type == allele_converted)]
                df_BA = add_column(df_BA, 'BA-score', data, records, gene, protein, length, allele, mutations, strain, sample)
                df_EL = add_column(df_EL, 'EL-score', data, records, gene, protein, length, allele, mutations, strain, sample)
        df_BA.columns = pd.MultiIndex.from_tuples(df_BA.columns, names=['sample', 'strain', 'mutation', 'allele'])
        df_BA.to_pickle(join('output/netMHCpan/pkl', protein + '_' + str(length) + 'mers_BA.pkl'))
        df_EL.columns = pd.MultiIndex.from_tuples(df_EL.columns, names=['sample', 'strain', 'mutation', 'allele'])
        df_EL.to_pickle(join('output/netMHCpan/pkl', protein + '_' + str(length) + 'mers_EL.pkl'))