from os import listdir, makedirs
from os.path import join, exists
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm

def add2list(peptide_dict, gene, tmp):
    peptides = peptide_dict[gene]
    for x in tmp:
        peptides.append(x)
    peptide_dict[gene] = list(pd.Series(peptides, dtype='object').unique())
    return peptide_dict

# peptide libraries
fasta = [f for f in listdir(join('input', 'fasta')) if 'fasta' in f]
for tq, length in zip(tqdm(range(8, 26), 'Preparing UM peptides'), range(8, 26)):
    peptide_dict = {}
    for f in fasta:
        protein = f.split('.')[0]
        peptide_dict[protein] = []
        records = [r for r in SeqIO.parse(join('input/fasta', f), 'fasta')]
        for r in records:
            tmp = [str(r.seq[i:i + length]) for i in range(len(r.seq) - length - 1) if 'X' not in str(r.seq[i:i+8])]
            peptide_dict = add2list(peptide_dict, protein, tmp)
    if exists(join('input', 'UM_peptides')) == False:
        makedirs(join('input', 'UM_peptides'))
    for protein in list(peptide_dict):
        peptides = peptide_dict[protein]
        with open(join('input/UM_peptides', protein + '_' + str(length) + 'mer.pep'), 'w') as output:
            for p in peptides:
                output.write(f"{p}\n")

# allele lists
print('Preparing MHCI allele list...')
df = pd.read_csv(join('input', 'data.csv'))

mhci_alleles = []
for c in ['A-1', 'A-2', 'B-1', 'B-2', 'C-1', 'C-2']:
    for a in df[c]:
        a = a.strip(' ').strip('g')
        out = 'HLA-' + a[0] + a.split('*')[1]
        if out not in mhci_alleles:
            mhci_alleles.append(out)
            
with open(join('input', 'MHCI_alleles.txt'), 'w') as output:
    for m in mhci_alleles:
        output.write(f"{m}\n")

print('Preparing MHCII allele list...')
mhcii_alleles = []
for c in ['DR-1', 'DR-2', 'DRw-1', 'DRw-2']:
    for a in df[c]:
        if isinstance(a, str):
            a = ':'.join(a.strip(' ').strip('g').split(' ')[0].split(':')[:2])
            out = a[:4] + '_' + a.split('*')[1].replace(':', '')
            if out not in mhcii_alleles:
                mhcii_alleles.append(out)
            
for c in ['DQB-1', 'DQB-2', 'DP-1', 'DP-2']:
    for i in ['1', '2']:
        if 'Q' in c:
            tmp = df[[c, c.replace('B', 'A')[:-1] + i]]
        elif 'P' in c:
            tmp = df[[c, c[:-2] + 'A-' + i]]
        for r in tmp.iterrows():
            b = r[1][0]
            b = ':'.join(b.strip(' ').strip('g').split(' ')[0].split(':')[:2]).replace(':', '').replace('*', '').split('(')[0].strip('Q')
            a = r[1][1]
            a = ':'.join(a.strip(' ').strip('g').split(' ')[0].split(':')[:2]).replace(':', '').replace('*', '').split('(')[0].strip('Q')
            out = 'HLA-' + a + '-' + b
            if out not in mhcii_alleles:
                mhcii_alleles.append(out)
            
with open(join('input', 'MHCII_alleles.txt'), 'w') as output:
    for m in mhcii_alleles:
        output.write(f"{m}\n")
        
print('Done!')