"""
Functions for identifying functional groups close to a specific substructure
within a molecule.  Functional groups are defined using Ertl's definition or 
using the list of predefined functional groups in RDKit.
"""


from rdkit import Chem
import numpy as np
import pandas as pd
import re

import sys, os

# Import functional groups code from RDKit Contrib:
# Add Contrib directory to path, this doesn't work on Archie so hardcoded 
# below instead:
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'IFG'))
#sys.path.insert(0, '/users/xpb20111/.conda/envs/deepchem/share/RDKit/Contrib/IFG')
from ifg import identify_functional_groups

from struct_fns import calc_bond_distance_between_atom_idxs, \
    get_substructure_atom_idxs_connected_to_mol
from brics_fns import cap_attachment_points, rm_attachment_point_charge, \
    get_attachment_point_numbers
from join_mols import join_fragments


#################################################
# Code using Ertl functional groups definition: #
#################################################


def get_nearest_fn_grps_to_core(cmpd, substruct):
    """
    Function to find nearest functional groups to specific substructure within
    a molecule, based on Ertl's definition of functional groups
    (https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0225-z).

    Parameters:
    cmpd (str or rdkit mol): Molecule
    substruct (str or rdkit mol): Substructure

    Returns:
    list of tuples: (Functional group structure, 
                     Distance to substructure (in bonds))

    >>> get_nearest_fn_grps_to_core(
    ...     'CCC(N)CC(c1ccccc1)CCOCCCCC(=O)OCCCCc1ccc(C)cc1',
    ...     'c1ccccc1')
    [('CN', 4), ('COC(C)=O', 5)]
    """

    if isinstance(cmpd, str):
        cmpd = Chem.MolFromSmiles(cmpd)
    if isinstance(substruct, str):
        substruct = Chem.MolFromSmiles(substruct)

    nearest_fn_grps = []

    fn_grps = identify_functional_groups(cmpd)

    substruct_atom_idxs = cmpd.GetSubstructMatches(substruct)

    fn_grp_types = [fn_grp.type for fn_grp in fn_grps]
    fn_grp_idxs = [fn_grp.atomIds for fn_grp in fn_grps]

    for substruct_atom_idxs in cmpd.GetSubstructMatches(substruct):

        # Still need to include connection points:
        for substruct_attachment_atom_idx in \
            get_substructure_atom_idxs_connected_to_mol(
                cmpd, substruct_match=substruct_atom_idxs):

            n_bonds = calc_bond_distance_between_atom_idxs(
                                            cmpd,
                                            [substruct_attachment_atom_idx],
                                            fn_grp_idxs,
                                            excluded_paths=substruct_atom_idxs)

            # If none of the functional groups have been found:
            if np.all([np.isnan(i) for i in n_bonds]):
                break

            nearest_fn_grp_i = np.argsort(n_bonds)[0]

            nearest_fn_grp = np.array(fn_grp_types)[nearest_fn_grp_i]
            dist_to_fn_grp = n_bonds[nearest_fn_grp_i]

            nearest_fn_grps.append((nearest_fn_grp, int(dist_to_fn_grp)))

    return nearest_fn_grps


############################################
# Code using predefined functional groups: #
############################################


def get_prefined_functional_groups():
    """
    Read functional groups from RDKit .txt file and make some modifications/
    additions to the definitions.  This function only needs to be run once to 
    generate the predefined_functional_groups.csv file.

    Generates file: predefined_functional_groups.csv
    """

    fn_grp_ls = []

    from rdkit.Data import __path__ as rdkit_data_path
    rdkit_data_path = rdkit_data_path._path[0]
    for line in open(rdkit_data_path+'/FunctionalGroups.txt', 'r').readlines():
    #for line in open('/users/xpb20111/.conda/envs/deepchem/share/RDKit/Data/FunctionalGroups.txt', 'r').readlines():

        line = line.strip()
        if line.startswith('//') or (line == ''): #startswith('\n'):
            continue
        else:
            fn_grp_info = line.split()

            if len(fn_grp_info) > 3:
                fn_grp_info[2] = '_'.join(fn_grp_info[2:])
                del fn_grp_info[3:]

            _, sma, fn_grp_name = fn_grp_info

            # Only include single bonds:
            if sma[1] != '-':
                continue

            try:
                fn_grp = Chem.MolToSmiles(Chem.MolFromSmarts(sma))
            except:
                print('Error: {}'.format(fn_grp_name))
                continue

            # Modify specific groups:
            if fn_grp_name == 'halogens':
                fn_grp_ls.append(('*F', 'fluorine'))
                fn_grp_ls.append(('*Cl', 'chlorine'))
                fn_grp_ls.append(('*Br', 'bromine'))
                fn_grp_ls.append(('*I', 'iodine'))

            elif fn_grp_name == 'nitro':
                fn_grp_ls.append(('*[N+]([O-])=O', 'nitro'))

            elif fn_grp_name == 'diazo':
                fn_grp_ls.append(('*N=N', 'diazo'))

            else:
                fn_grp_ls.append((fn_grp, fn_grp_name))

    # Add additional functional groups:
    fn_grp_ls.insert(1, ('*NC(N)=O', 'urea'))
    
    with open('predefined_functional_groups.csv', 'w') as out:
        out.write('SMARTS,Name\n')
        out.write('\n'.join([','.join(i) for i in fn_grp_ls]))


def get_predefined_fn_grp_core_frgs(cmpd, 
                                    substruct, 
                                    fn_grp_ls='predefined_functional_groups.csv'):
    """
    Function to find predefined functional groups immediately adjacent to 
    specific substructures within a molecule.

    Parameters:
    cmpd (str or rdkit mol): Molecule
    substruct (str or rdkit mol): Substructure
    fn_grp_ls: List of functional groups from reading 
               predefined_functional_groups.csv

    Returns:
    List of tuples: (Functional group structure,
                     Functional group bonded to substructure,
                     Functional group name)
    
    >>> smiles = 'CC(O)(CS(=O)(=O)c1ccc(F)cc1)C(=O)Nc1ccc(C#N)c(C(F)(F)F)c1'
    >>> substruct = 'FC(F)(F)c1ccccc1'
    >>> get_predefined_fn_grp_core_frgs(smiles, substruct, fn_grp_ls='predefined_functional_groups.csv') #, canonicalise=True)
    [['*C#N', 'N#Cc1ccccc1C(F)(F)F', 'cyano'], ['*NC(C)=O', 'CC(=O)Nc1cccc(C(F)(F)F)c1', 'methyl_amide']]
    
    >>> smiles = 'O=C(NNc1c([N+](=O)[O-])cc(C(F)(F)F)cc1[N+](=O)[O-])c1ccco1'
    >>> get_predefined_fn_grp_core_frgs(smiles, substruct, fn_grp_ls='predefined_functional_groups.csv')
    [['*[N+]([O-])=O', '[H]c1c([H])cc(C(F)(F)F)cc1[N+](=O)[O-]', 'nitro'], \
['*[N+]([O-])=O', '[H]c1cc(C(F)(F)F)cc([N+](=O)[O-])c1[H]', 'nitro'], \
['*N', '[H]c1cc(C(F)(F)F)cc([H])c1N', 'primary_amines']]
    """

    fn_grps = []

    if isinstance(cmpd, str):
        cmpd = Chem.MolFromSmiles(cmpd)
    if isinstance(substruct, str):
        substruct = Chem.MolFromSmiles(substruct)

    if isinstance(fn_grp_ls, str):
        fn_grp_ls = pd.read_csv(fn_grp_ls).to_numpy()

    core_with_attachment = Chem.ReplaceSidechains(cmpd, substruct)

    #Chem.Kekulize(core_with_attachment)
    core_with_attachment_smi = Chem.MolToSmiles(core_with_attachment) #, kekuleSmiles=True)
    
    # Check that core is a valid SMILES string:
    if Chem.MolFromSmiles(core_with_attachment_smi) is None:
        return [(None, None, None)]

    # Remove charge associated with attachment points:
    core_with_attachment_smi = rm_attachment_point_charge(core_with_attachment_smi)

    for at_point in get_attachment_point_numbers(core_with_attachment_smi):
        for fn_grp, fn_grp_name in fn_grp_ls:
            fn_grp_at = re.sub('\*', '['+str(at_point)+'*]',
                            fn_grp)

            substruct_frg = join_fragments(fn_grp_at, core_with_attachment_smi)
            substruct_frg = cap_attachment_points(substruct_frg, canonicalise=True) #False, check_valid=False)

            substruct_frg_mol = Chem.MolFromSmiles(substruct_frg)

            if substruct_frg_mol is not None:

                if cmpd.HasSubstructMatch(substruct_frg_mol):
                    fn_grps.append((fn_grp, substruct_frg, fn_grp_name))
                    # Break here due to functional group hierarchy?
                    break

    return fn_grps
