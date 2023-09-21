from rdkit import Chem
#from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina
import numpy as np


def get_butina_clustering(smiles=[], fp=[], cutoff=0.6):
    """
    Do Butina clustering.
    """
    if len(fp) == 0:
        fp = [Chem.RDKFingerprint(Chem.MolFromSmiles(smi))        
              for smi in smiles]

    dist = []
    for i in range(1, len(fp)):
        # Lower triangle:
        sim = DataStructs.BulkTanimotoSimilarity(fp[i],fp[:i])
        dist.extend([1 - s for s in sim])

    clusters = Butina.ClusterData(dist, len(fp), cutoff, isDistData=True)

    cluster_membership = np.zeros(len(fp), dtype=int)
    for cl_i, cl in enumerate(clusters):
        cluster_membership[list(cl)] = cl_i
    return cluster_membership
