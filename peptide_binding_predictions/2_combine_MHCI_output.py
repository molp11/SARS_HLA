from os import listdir
from os.path import join
import pandas as pd
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

files = [f for f in listdir('output/netMHCpan') if f.split('.')[-1] == 'xls']
master = pd.DataFrame()
for tq, file in zip(tqdm(range(len(files)), 'Combining netMHCpan output'), files):
    df = pd.read_table(join('output/netMHCpan', file), header=1, delim_whitespace=True)
    df['Protein'] = [file.split('_')[0]]*len(df)
    df['Peptide_length'] = [file.split('_')[1].strip('mer')]*len(df)
    df['HLA_type'] = [file.split('_')[2].strip('.xls')]*len(df)
    df['Peptide_length'] = [file.split('_')[1].strip('mer')]*len(df)
    if len(master) == 0:
        master = df
    else:
        master = pd.concat([master, df])
        
master.to_csv(join('output/netMHCpan', 'NetMHCpan-output.csv'), index=False)