"""
Functions to join molecular fragments with a single bond.
"""


from rdkit import Chem
import re


def join_fragments(smi_1, smi_2):
    mol = None
    for join_point in re.findall('\[[0-9]\*\]', smi_1):
        if join_point in smi_2:
            mol = Chem.MolFromSmiles(smi_1+'.'+smi_2)
            if mol is not None:
                mol = make_bond(mol, join_idx=int(join_point[1]))
            break
    if mol is None:
        print('ERROR')
    return mol


def make_bond(mol, join_idx=None, returnMol=False):
    if join_idx is not None:
        join_idx = int(join_idx)

    rm_atoms = []
    join_atoms = []
    for i, at in enumerate(mol.GetAtoms()):
        if (at.GetSymbol() == '*') and ((join_idx is None) or (at.GetIsotope() == join_idx)):
            rm_atoms.append(at.GetIdx())
            join_atoms.append(at.GetNeighbors()[0].GetIdx())

    # Maybe use this: https://sourceforge.net/p/rdkit/mailman/message/31703978/
    # Editable mol, see: https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html
    ed_mol = Chem.EditableMol(mol)
    # Need to do in batch:
    ed_mol.BeginBatchEdit()
    for at in rm_atoms:
        ed_mol.RemoveAtom(at)
    ed_mol.AddBond(*join_atoms, Chem.BondType.SINGLE)
    ed_mol.CommitBatchEdit()
    # Convert back from editable mol:
    mol = ed_mol.GetMol()

    Chem.SanitizeMol(mol)

    if returnMol:
        return mol
    else:
        return Chem.MolToSmiles(mol)
