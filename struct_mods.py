"""
Functions to modify molecular structures.
"""


from rdkit import Chem
import sys
sys.path.insert(0, '/users/xpb20111/programs')
from cheminfo_utils.brics_fns import get_attachment_point_numbers
from cheminfo_utils.join_mols import make_bond
import re


def swap_side_chains(cmpd, 
                     core, 
                     calc_mass_imbalance=False, 
                     calc_heavy_atom_imbalance=False):
    """
    Swap positions of two side chains around a core.
    
    >>> new_mol, side_chain_imbalance = swap_side_chains(
    ...     cmpd='CCCCCCN1CN(c2cccc3c[nH]nc23)C(c2ccccc2)C1=O', 
    ...     core='N1CNC(c2ccccc2)C1=O', 
    ...     calc_heavy_atom_imbalance=True)
    >>> Chem.MolToSmiles(new_mol), side_chain_imbalance
    ('CCCCCCN1CN(c2cccc3c[nH]nc23)C(=O)C1c1ccccc1', 3)
    """

    #if isinstance(cmpd, str):
    #    cmpd = Chem.MolFromSmiles(cmpd)
    #
    #if isinstance(core, str):
    #    core = Chem.MolFromSmiles(core)
 
    sidechains = Chem.ReplaceCore(cmpd, core)
    if sidechains is None:
        print('WARNING: Cannot find core in molecule: {}'.format(
            Chem.MolToSmiles(cmpd)))
        if calc_heavy_atom_imbalance:
            return None, None
        else:
            return None


    sidechains = Chem.MolToSmiles(sidechains)

    # Check there are two side chains:
    if '.' not in sidechains:
        print('WARNING: Only one side chain found, cannot swap side chains.')
        if calc_heavy_atom_imbalance:
            return None, None
        else:
            return None

    if calc_heavy_atom_imbalance:
        heavy_atoms = ['c', 'n', 'o', 's', 'f', 'cl', 'br', 'i']

        separate_chains = []
        for sc in sidechains.split('.'):
            at_point = get_attachment_point_numbers(sc)
            if len(at_point) != 1:
                raise ValueError()
            separate_chains.append([at_point[0], 
                                    len([atom for atom in sc.lower() 
                                         if atom in heavy_atoms])])

        separate_chains = sorted(separate_chains, key=lambda x: x[0])
        heavy_atom_imbalance = separate_chains[1][1] - separate_chains[0][1]

    if calc_mass_imbalance:
        raise NotImplmenetedError()

    core = Chem.MolToSmiles(Chem.ReplaceSidechains(cmpd, core))
    
    at_points = get_attachment_point_numbers(sidechains)
    if len(at_points) != 2:
        raise ValueError(at_points)
    
    # Swap attachment point numbering:
    sidechains = re.sub('\['+str(at_points[0])+'\*\]', '[_*]', sidechains)
    sidechains = re.sub('\['+str(at_points[1])+'\*\]', '['+str(at_points[0])+'*]', sidechains)
    sidechains = re.sub('\[\_\*\]', '['+str(at_points[1])+'*]', sidechains)
    
    mol = Chem.MolFromSmiles(core+'.'+sidechains)
    for at_n in at_points:
        mol = make_bond(mol, int(at_n), returnMol=True)

    if calc_heavy_atom_imbalance:
        return mol, heavy_atom_imbalance
    else:
        return mol
