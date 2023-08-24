from rdkit import Chem
import numpy as np


# Taken from: ~/pIC50/notebooks/Process_PDB_ligands_from_Ellen.ipynb
def GetSubstructRMSD(mol1, mol2, substruct):
    """
    Use spyrmsd to get substructure RMSD 
    taking symmetry into account.
    """

    from spyrmsd import io, rmsd, molecule

    mols = [mol1, mol2]
    spyrmsd_mol = [None, None]
    coords_spyrmsd = [None, None]
    anum_spyrmsd = [None, None]
    adj_spyrmsd = [None, None]
    atomIds = [None, None]
    coords_spyrmsd_sub = [None, None]
    anum_spyrmsd_sub = [None, None]
    adj_spyrmsd_sub = [None, None]

    best_substruct_rmsd = 1000000
    
    for mol_i, mol in enumerate(mols):
        spyrmsd_mol[mol_i] = molecule.Molecule.from_rdkit(mol)
        spyrmsd_mol[mol_i].strip()
        coords_spyrmsd[mol_i] = spyrmsd_mol[mol_i].coordinates
        anum_spyrmsd[mol_i] = spyrmsd_mol[mol_i].atomicnums
        adj_spyrmsd[mol_i] = spyrmsd_mol[mol_i].adjacency_matrix

        atomIds[mol_i] = mol.GetSubstructMatches(substruct)

        if len(atomIds[mol_i]) == 0:
            print('WARNING: No substructure match found in compound: {}'.format(mol_i))
            best_substruct_rmsd = np.nan
            return best_substruct_rmsd

        elif len(atomIds[mol_i]) > 1:
            print('WARNING: Compound {} has {} substructure matches, will return smallest RMSD'.format(mol_i, len(atomIds_mol1)))

    for atomIds_sub1 in atomIds[0]:
        atomIds_sub1 = list(atomIds_sub1)
        coords_spyrmsd_sub[0] = coords_spyrmsd[0][atomIds_sub1]
        anum_spyrmsd_sub[0] = anum_spyrmsd[0][atomIds_sub1]
        adj_spyrmsd_sub[0] = adj_spyrmsd[0][atomIds_sub1][:,atomIds_sub1]

        for atomIds_sub2 in atomIds[1]:
            atomIds_sub2 = list(atomIds_sub2)
            coords_spyrmsd_sub[1] = coords_spyrmsd[1][atomIds_sub2]
            anum_spyrmsd_sub[1] = anum_spyrmsd[1][atomIds_sub2]
            adj_spyrmsd_sub[1] = adj_spyrmsd[1][atomIds_sub2][:,atomIds_sub2]

            substruct_rmsd_sym = rmsd.symmrmsd(*coords_spyrmsd_sub,
                                               *anum_spyrmsd_sub,
                                               *adj_spyrmsd_sub)

            if substruct_rmsd_sym < best_substruct_rmsd:
                best_substruct_rmsd = substruct_rmsd_sym
    
    return best_substruct_rmsd


