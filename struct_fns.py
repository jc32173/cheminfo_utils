"""
Functions for calculating structural features of molecules.
"""


from rdkit import Chem
import numpy as np
import collections


def calc_bond_distance_between_substructs(cmpd, 
                                          sub1, 
                                          sub2):
    """
    Calculate the distance in bonds between substructures in a molecule.

    Parameters:
    cmpd (str or rdkit mol): Molecule
    sub1 (str or rdkit mol): First substructure
    sub2 (list of str or rdkit mol): List of second substructures

    Returns:
    list: Minimum distance in bonds from first substructure to each 
          second substructure in molecule

    >>> calc_bond_distance_between_substructs('CCC(N)CCCCOCCCCC(=O)OCCCCc1ccccc1',
    ...                                       'c1ccccc1',
    ...                                       ['CNCN', 'COC', 'CN', 'COC(C)=O'])
    [nan, 4.0, 16.0, 4.0]
    """

    # Assume any string inputs are SMILES and convert to
    # RDKit mol objects:
    if isinstance(cmpd, str):
        cmpd = Chem.MolFromSmiles(cmpd)
    if isinstance(sub1, str):
        sub1 = Chem.MolFromSmiles(sub1)
    if isinstance(sub2, str):
        sub2 = [sub2]
    for sub_i, sub in enumerate(sub2):
        if isinstance(sub, str):
            sub2[sub_i] = Chem.MolFromSmiles(sub)

    # Get atom indices for first substructure:
    atom_ring = cmpd.GetSubstructMatches(sub1)
    if len(atom_ring) == 0:
        raise ValueError('Substructure 1 not found in molecule ({})'.format(Chem.MolToSmiles(sub1)))
    atom_ring = tuple([i for j in atom_ring for i in j])

    # Get atom indices for second substructures:
    sub2_idxs = []
    for sub_i, sub in enumerate(sub2):
        sub2_idx = cmpd.GetSubstructMatches(sub)
        if len(sub2_idx) == 0:
            #print('WARNING: Substructure 2 not found in molecule ({})'.format(Chem.MolToSmiles(sub)))
            sub2_idxs.append((np.nan,)) #-1
            sub2[sub_i] = None
            #found_dist[sub_i] = True
            #n_bonds[sub_i] = -1
        else:
            sub2_idx = tuple([i for g in sub2_idx for i in g])
            sub2_idxs.append(sub2_idx)

    return calc_bond_distance_between_atom_idxs(cmpd, atom_ring, sub2_idxs, substructs=sub2)


def calc_bond_distance_between_atom_idxs(cmpd, 
                                         atoms1, 
                                         atoms2, 
                                         substructs=None, 
                                         excluded_paths=[]):
    """
    Calculate the shortest distance in bonds between groups of atom indices 
    in a molecule.

    Parameters:
    cmpd (str or rdkit mol): Molecule
    atoms1 (str or rdkit mol): First group of atom indices
    atoms2 (list of str or rdkit mol): List of second groups of atom indices
    substructs (list of str): IDs or SMILES for groups of atoms in sub2_idx, only used for printing progress

    Returns:
    list: Minimum distance in bonds from first group of atom indices to each second group in the molecule

    >>> calc_bond_distance_between_atom_idxs('CCC(N)CCCCOCCCCC(=O)OCCCCc1ccccc1',
    ...                                      (20, 21, 22, 23, 24, 25),
    ...                                      [(26, 27), (3, ), (8, ), (13, 14, 15)])
    [nan, 17.0, 11.0, 5.0]
    """

    # Assume string input is SMILES and convert to
    # RDKit mol object:
    if isinstance(cmpd, str):
        cmpd = Chem.MolFromSmiles(cmpd)

    # Ensure atoms2 has dimension == 2, since it should be a list of 
    # lists of atom indices:
    if not isinstance(atoms2, collections.abc.Iterable):
        atoms2 = [atoms2]
    if not isinstance(atoms2[0], collections.abc.Iterable):
        atoms2 = [atoms2]

    atom_ring = atoms1

    checked_atoms = []
    found_dist = np.array([False]*len(atoms2))
    n_bonds = np.array([0]*len(atoms2), dtype=float)

    # Expand from first group of atoms along bonds until all other
    # atom groups are found:
    while True:
        for atom_grp_i, atom_idxs in enumerate(atoms2):
            if not found_dist[atom_grp_i]:
                if len(set(atom_idxs) & set(atom_ring)) > 0:
#                    if substructs is None:
#                        print('Found substructure: {} after {} bonds'.format(
#                              substruct[atom_grp_i], n_bonds[atom_grp_i]))
#                    else:
#                        print('Found atom indices: {} after {} bonds'.format(
#                              atom_idx, n_bonds[atom_grp_i]))
                    found_dist[atom_grp_i] = True

        # Stop if all substructures have been found:
        if np.all(found_dist):
            break

        checked_atoms += atom_ring
        prev_atom_ring = atom_ring
        atom_ring = []
        n_bonds[~found_dist] += 1

        # Otherwise move out from first substucture by 1 bond:
        for at in prev_atom_ring:
            nns = cmpd.GetAtomWithIdx(at).GetNeighbors()
            for nn in nns:
                nn_idx = nn.GetIdx()
                if (nn_idx not in checked_atoms) and \
                   (nn_idx not in excluded_paths):
                    atom_ring.append(nn_idx)
        
        # Stop if all atoms have been checked:
        if len(atom_ring) == 0:
            n_bonds[~found_dist] = np.nan #-1
            break

    return list(n_bonds)


def calc_path_between_atom_idxs(cmpd, atoms1, atoms2):
    """
    Function to calculate the shortest path through a molecule which connects 
    specific atoms.
    """
    raise NotImplementedError()
#    # Assume string input is SMILES and convert to
#    # RDKit mol object:
#    if isinstance(cmpd, str):
#        cmpd = Chem.MolFromSmiles(cmpd)
#
##    # Ensure atoms2 has dimension == 2, since it should be a list of 
##    # lists of atom indices:
##    if not isinstance(atoms2, collections.abc.Iterable):
##        atoms2 = [atoms2]
##    if not isinstance(atoms2[0], collections.abc.Iterable):
##        atoms2 = [atoms2]
#
#    atom_ring = atoms1
#
#    checked_atoms = []
#    found_dist = np.array([False]*len(atoms2))
#    n_bonds = np.array([0]*len(atoms2), dtype=float)
#
#    # Expand from first group of atoms along bonds until all other
#    # atom groups are found:
#    while True:
#        for atom_grp_i, atom_idxs in enumerate(atoms2):
#            if not found_dist[atom_grp_i]:
#                if len(set(atom_idxs) & set(atom_ring)) > 0:
##                    if substructs is None:
##                        print('Found substructure: {} after {} bonds'.format(
##                              substruct[atom_grp_i], n_bonds[atom_grp_i]))
##                    else:
##                        print('Found atom indices: {} after {} bonds'.format(
##                              atom_idx, n_bonds[atom_grp_i]))
#                    found_dist[atom_grp_i] = True
#
#        # Stop if all substructures have been found:
#        if np.all(found_dist):
#            break
#
#        checked_atoms += atom_ring
#        prev_atom_ring = atom_ring
#        atom_ring = []
#        n_bonds[~found_dist] += 1
#
#        # Otherwise move out from first substucture by 1 bond:
#        for at in prev_atom_ring:
#            nns = cmpd.GetAtomWithIdx(at).GetNeighbors()
#            for nn in nns:
#                nn_idx = nn.GetIdx()
#                if (nn_idx not in checked_atoms) and \
#                   (nn_idx not in excluded_paths):
#                    atom_ring.append(nn_idx)
#        
#        # Stop if all atoms have been checked:
#        if len(atom_ring) == 0:
#            n_bonds[~found_dist] = np.nan #-1
#            break
#
#    return list(n_bonds)


def get_substructure_atom_idxs_connected_to_mol(cmpd, substruct=None, substruct_match=None):
    """
    Get the indices of the atoms in a substructure which are connected to the 
    rest of the molecule.

    Parameters:
    cmpd (str or rdkit mol): Molecule
    substruct (str or rdkit mol, opt.): Substructure
    substruct_match (list or tuple, opt.): Atom indices of substructure match

    Returns:
    list: Atom indices of atoms in the substructure
          which are connected to the rest of the molecule

    >>> get_substructure_atom_idxs_connected_to_mol(
    ...     'CCC(N)CC(c1ccccc1)CCOCCCCC(=O)OCCCCc1ccc(C)cc1',
    ...     'c1ccccc1')
    [6, 26, 29]
    """

    if (substruct is None) and (substruct_match is None):
        raise ValueError

    atom_idxs = []

    if isinstance(cmpd, str):
        cmpd = Chem.MolFromSmiles(cmpd)

    if substruct_match is not None:
        substruct_atom_matches = [substruct_match]
    elif isinstance(substruct, str):
        substruct = Chem.MolFromSmiles(substruct)
    if substruct_match is None:
        substruct_atom_matches = cmpd.GetSubstructMatches(substruct)

    for substruct_atom_idxs in substruct_atom_matches:
        substruct_atom_idxs = set(substruct_atom_idxs)
        for atom_idx in substruct_atom_idxs:
            nn_idxs = [nn.GetIdx() for nn in cmpd.GetAtomWithIdx(atom_idx).GetNeighbors()]
            if len(set(nn_idxs) - substruct_atom_idxs) > 0:
                atom_idxs.append(atom_idx)

    return atom_idxs
