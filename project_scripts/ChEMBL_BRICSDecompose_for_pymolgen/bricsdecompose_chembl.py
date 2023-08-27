import multiprocessing as mp
import numpy as np
import pandas as pd
from datetime import datetime
import gzip
import sys, os
from rdkit import Chem
from rdkit.Chem.BRICS import BRICSDecompose
from rdkit.Chem.SaltRemover import SaltRemover

#sys.path.insert(-1, '/users/xpb20111/scripts/sys/')
#from slurm.slurm_stats import slurm_stats
#get_slurm_stats = slurm_stats(filename='slurm_stats.csv', flush_file=True)


def run_bricsdecompose(smi, chembl_id=None):
    start_t = datetime.now()
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print('ERROR: Cannot read: {} {}'.format(chembl_id, smi))
        brics_frgms = []
    else:
        mol_stripped = saltremover.StripMol(mol, dontRemoveEverything=True)
        brics_frgms = list(BRICSDecompose(mol_stripped, 
                                          keepNonLeafNodes=False, 
                                          singlePass=False))
    end_t = datetime.now()
    time_taken = (end_t - start_t).total_seconds()
    return brics_frgms, time_taken


# Get start_row and n_rows from command line:
start_line = int(sys.argv[1])
end_line = int(sys.argv[2])
#array_job_id=int(sys.argv[1])
array_job_id = os.environ['SLURM_ARRAY_TASK_ID']

n_lines = end_line - start_line
if n_lines == 0:
    raise ValueError('start_line == end_line, will not read any lines from file')

chembl_file = "/users/xpb20111/ChEMBL/chembl_33_chemreps.txt.gz"

#brics_outfile = 'chembl_33_brics_{}.txt.gz'.format(int(start_line/n_lines))
#frgs_outfile = 'chembl_33_brics_frgs_{}.txt'.format(int(start_line/n_lines))
brics_outfile = 'chembl_33_brics_{}.txt.gz'.format(array_job_id)
frgs_outfile = 'chembl_33_brics_frgs_{}.txt'.format(array_job_id)
frgs_out = open(frgs_outfile, 'a')

saltremover = SaltRemover()

# Get number of molecules from number of lines in chembl file:
with gzip.open(chembl_file, 'rb') as f:
    n_mols = sum(1 for _ in f) - 1

# Iterate 1 line at a time:
df_iter = pd.read_csv(chembl_file, 
                      header=0, 
                      # Skip lines after header:
                      skiprows=range(1, start_line), 
                      nrows=n_lines, 
                      chunksize=1, 
                      sep='\t', 
                      usecols=['chembl_id', 'canonical_smiles'])

print('Running BRICSDecompose...')

write_header = True
if os.path.isfile(brics_outfile):
    write_header = False

for f_i, df in enumerate(df_iter):
    #get_slurm_stats()
    smiles = df['canonical_smiles'].squeeze()
    brics_frgms, time_taken = run_bricsdecompose(smiles)
    df['brics_fragments'] = [brics_frgms]
    df['bricsdecompose_time'] = time_taken
    df.to_csv(brics_outfile, 
              sep='\t', mode='a', header=write_header, index=False)
    write_header = False

    if frgs_outfile != '':
        canon_frgs = [Chem.MolToSmiles(Chem.MolFromSmiles(f)) 
                      for f in df['brics_fragments'].squeeze()]
        frgs_out.write('\n'.join(canon_frgs)+'\n')
        frgs_out.flush()
