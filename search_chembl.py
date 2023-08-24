from chembl_webresource_client.new_client import new_client
import pandas as pd


# Function to search ChEMBL for compounds with values for 
# particular properties and with similarity to query compound.
# Taken from notebook: 
def search_chembl(smiles, val_type, similarity_cutoff):
    """
    Function to search ChEMBL
    """
    if similarity_cutoff < 40:
        print('ERROR: Minimum similarity score: 40% (see error message from python)')
        return None
    # Do search for molecules within given similarity cutoff:
    similarity = new_client.similarity
    mols = similarity.filter(smiles=smiles, similarity=similarity_cutoff)\
                     .only(['molecule_chembl_id',
                            'similarity',
                            'molecule_structures'])
    print('Number of molecules found: {} (similarity: {}%)'.format(len(mols), similarity_cutoff))
    # Do search for activities of type "LogD" within subset of molecules found above:
    activity = new_client.activity
    activities = activity.filter(molecule_chembl_id__in=[m['molecule_chembl_id'] for m in mols])\
                         .filter(standard_type=val_type)\
                         .only(['molecule_chembl_id',
                                'target_chembl_id',
                                'value',
                                'type',
                                'relation',
                                'units',
                                'assay_description'])
    if len(activities) == 0:
        print('No activities found for these molecules for {}.'.format(val_type))
        return None
    # Combine previous two searches into a single DataFrame:
    df = pd.merge(left=pd.DataFrame(mols), left_on='molecule_chembl_id', 
                  right=pd.DataFrame(activities), right_on='molecule_chembl_id',
                  how='outer')
    df['canonical_smiles'] = df['molecule_structures'].apply(lambda row: row['canonical_smiles'])
    df.drop(columns='molecule_structures', inplace=True)
    return df
