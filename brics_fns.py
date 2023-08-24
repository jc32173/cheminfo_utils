from rdkit import Chem
from rdkit.Chem.BRICS import BRICSDecompose
import numpy as np
import re
#from rdkit_funcs import canonicalise_tautomer
#from smi_funcs import canonicalise_tautomer


# Don't really need canonicalise_frgms and rm_BRICS_rules, but
# may be needed for old code:

def canonicalise_frgms(smi, rm_BRICS_rules=False):
    """
    Canonicalise BRICS fragment.
    """
    
    if rm_BRICS_rules:
        smi = re.sub('\[[0-9]+\*\]', '*', smi)
    smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    return smi


def rm_BRICS_rules(smi, canonicalise=True):
    """
    Remove BRICS numbers associated with attachment points.
    """

    smi = re.sub('\[[0-9]+\*\]', '*', smi)
    if canonicalise:
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    return smi


def rm_attachment_point_charge(smi):
    """
    Sometimes attachment points still have charge, e.g.:

    >>> smiles = 'O=C(NNc1c([N+](=O)[O-])cc(C(F)(F)F)cc1[N+](=O)[O-])c1ccco1'
    >>> substruct = 'FC(F)(F)c1ccccc1'
    >>> core = Chem.MolToSmiles(
    ...            Chem.ReplaceSidechains(Chem.MolFromSmiles(smiles), 
    ...                                   Chem.MolFromSmiles(substruct)))
    >>> core
    '[1*+]c1cc(C(F)(F)F)cc([3*+])c1[2*]'
    >>> rm_attachment_point_charge(core)
    '[1*]c1cc(C(F)(F)F)cc([3*])c1[2*]'
    """

    smi = re.sub('\*[1-9]*[\+-]', '*', smi)
    return smi


def cap_attachment_points(smi, canonicalise=True, check_valid=True):
    """
    Cap attachment points (labelled with *) with H.
    """

    smi = re.sub('\[[0-9]+\*\]', '*', smi)
    smi = re.sub('\(\*\)', '([H])', smi)
    smi = re.sub('\*', '[H]', smi)

    # Check RDKit can read modified SMILES
    # and canonicalise:
    if check_valid or canonicalise:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print('ERROR: RDKit cannot read capped SMILES')
            return None
        if canonicalise:
            return Chem.MolToSmiles(mol)
        else:
            return smi
    else:
        return smi


def get_attachment_point_numbers(smi):
    """
    Get list of attachment point numbers.
    """

    at_point_nos = re.findall('\[([0-9]+)\*\]', smi)
    at_point_nos = [int(n) for n in at_point_nos]
    return at_point_nos


def get_atom_idxs_with_attachment_points(cmpd, attachment_point_numbers=[]):
    """
    Get indices of atoms with attachment points.

    Parameters:
    cmpd (str or rdkit mol): Molecule
    attachment_point_numbers (, opt.): Specific attachment point numbers

    Returns:
    list: Atom indices of atoms with attachment points

    >>> get_atom_idxs_with_attachment_points('[12*]CCCC(=O)OCC[3*]', [3, 12, 13])
    [8, 1, None]
    """

    # Number of attachment points to search for:
    n_points = len(attachment_point_numbers)
    if n_points > 0:
        atom_idxs = [None]*n_points
    else:
        atom_idxs = []

    if isinstance(cmpd, str):
        cmpd = Chem.MolFromSmiles(cmpd)

    for at in cmpd.GetAtoms():
        for nn in at.GetNeighbors():
            if nn.GetSymbol() == '*':
                # Keep order of attachment points if specific numbers given:
                if n_points > 0:
                    # Attachment point number is set as isotope:
                    at_point_no = nn.GetIsotope()
                    if at_point_no in attachment_point_numbers:
                        atom_idxs[attachment_point_numbers.index(at_point_no)] = at.GetIdx()
                else:
                    atom_idxs.append(at.GetIdx())
    return atom_idxs


def rm_chirality(smi):
    """
    Remove chirality definition from a SMILES using regex.
    """

    smi = re.sub('\[C@+H?\]', 'C', smi)
    return smi


def core_BRICSDecompose(smi,
                        core='Cc1noc(C)c1c1ccccc1',
                        rm_BRICS_rules=False,
                        save_leaf_frags=False,
                        canonicalise_tautomers=False,
                        isomericSmiles=False):
    """
    Do BRICSDecompose to generate all possible intermediate fragments and organise the results
    based on whether they contain the specified core structure or are leaf nodes.
    """

    if isinstance(core, str):
        core = Chem.MolFromSmiles(core)

    mol = Chem.MolFromSmiles(smi)

    decomposed_smis = BRICSDecompose(mol, keepNonLeafNodes=True, singlePass=False)

    if canonicalise_tautomers:
        decomposed_smis = [canonicalise_tautomer(s) for s in decomposed_smis]

    if not isomericSmiles:
        decomposed_smis = [Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=False) for s in decomposed_smis]

    incomplete_frgs_core = []
    incomplete_frgs_no_core = []

    # for s in sorted(b_decom, key=len):
    # for s in sorted(b_decom, key=lambda x: len(re.sub('\[[0-9]\*\]', '', x))):

    # Order fragments by number of C, N, O atoms, determined by count in string:
    for frg_smi in sorted(decomposed_smis, key=lambda s: sum([s.lower().count(c) for c in ['c', 'n', 'o']])):
        if Chem.AddHs(Chem.MolFromSmiles(frg_smi)).HasSubstructMatch(core):
            incomplete_frgs_core.append(canonicalise_frgms(frg_smi, rm_BRICS_rules=rm_BRICS_rules))
        else:
            incomplete_frgs_no_core.append(canonicalise_frgms(frg_smi, rm_BRICS_rules=rm_BRICS_rules))

    if save_leaf_frags:
        leaf_smis = BRICSDecompose(mol, keepNonLeafNodes=False) #, singlePass=False)
        leaf_smis = [canonicalise_frgms(frg_smi, rm_BRICS_rules=rm_BRICS_rules) for frg_smi in leaf_smis]

        if canonicalise_tautomers:
            leaf_smis = [canonicalise_tautomer(s) for s in leaf_smis]

        if not isomericSmiles:
            leaf_smis = [Chem.MolToSmiles(Chem.MolFromSmiles(s), isomericSmiles=False) for s in leaf_smis]
    else:
        leaf_smis = []

    return incomplete_frgs_core, incomplete_frgs_no_core, leaf_smis


# Originally in: /users/xpb20111/notebooks/202210_Decomposed_cores_for_Bruno/202210_Partial_BRICS_decomposition.ipynb
def frgs_core_to_mol(smi, core_smi):
    """
    Get number of and SMILES of BRICS fragments between full molecule and core
    Note that the order of fragments is not preserved, use 
    get_common_core_intermediates if that is important.
    """

    # Original molecule:
    mol = Chem.MolFromSmiles(smi)

    # Substitute attachment points with H in core fragment:
    core_smi = re.sub('\(\*\)', '([H])', core_smi)
    core = Chem.MolFromSmiles(re.sub('\*', '[H]', core_smi))

    # Get side chains (difference between full molecule and core)
    # by removing core from molecule:
    side_chain_mol = Chem.ReplaceCore(mol, coreQuery=core)
    side_chain_smis = Chem.MolToSmiles(side_chain_mol).split('.')

    # Use BRICS to decompose side chains:
    n_frgs = 0
    frg_smis = []
    for s in side_chain_smis:
        n_frgs += len(BRICSDecompose(Chem.MolFromSmiles(s)))
        frg_smis += list(BRICSDecompose(Chem.MolFromSmiles(s)))

    return n_frgs, frg_smis


# Temporarily added function here, eventually import from smi_funcs.py:
def canonicalise_tautomer(cmpd):
    # Convert SMILES to RDKit Mol object:
    if isinstance(cmpd, str):
        cmpd = Chem.MolFromSmiles(cmpd)
    enumerator = rdMolStandardize.TautomerEnumerator()
    # enumerator.SetReassignStereo = True
    # enumerator.SetRemoveBondStereo = False
    # enumerator.SetRemoveSp3Stereo = False
    mol = enumerator.Canonicalize(cmpd)
    return mol


def get_common_core_intermediates(cmpd,
                                  core,
                                  only_smallest_n_fragments=None,
                                  cap_fragments=False,
                                  BRICS_only=True,
                                  rm_BRICS_rules=False,
                                  canonicalise_tautomers=False,
                                  isomericSmiles=False):
    """
    Get N BRICS fragments from a molecule which contain a 
    given core structure.
    """

    if isinstance(cmpd, str):
        mol = Chem.MolFromSmiles(cmpd)
    else:
        mol = cmpd

    if isinstance(core, str):
        core = Chem.MolFromSmiles(core)

    decomposed_mols = BRICSDecompose(mol,
                                     keepNonLeafNodes=True,
                                     singlePass=False,
                                     returnMols=True)
    decomposed_mols = np.array(list(decomposed_mols))

    incomplete_frgs_core = []

    # Order fragments by number of heavy atoms:
    decomposed_mols_n_atoms = np.array([m.GetNumHeavyAtoms() for m in decomposed_mols])

    n_atoms_core = core.GetNumHeavyAtoms()

    decomposed_mols = decomposed_mols[np.where(decomposed_mols_n_atoms >= n_atoms_core)]
    decomposed_mols_n_atoms = decomposed_mols_n_atoms[np.where(decomposed_mols_n_atoms >= n_atoms_core)]

    decomposed_mols = decomposed_mols[np.argsort(decomposed_mols_n_atoms)]
    decomposed_mols_n_atoms = np.sort(decomposed_mols_n_atoms)

    core_with_attachment = None
    for frg, n_atoms in zip(decomposed_mols, decomposed_mols_n_atoms):
        if Chem.AddHs(frg).HasSubstructMatch(core):
            if canonicalise_tautomers:
                frg = canonicalise_tautomer(frg)
            if not isomericSmiles:
                Chem.rdmolops.RemoveStereochemistry(frg)
            if cap_fragments:
                frg = Chem.MolFromSmiles(cap_attachment_points(Chem.MolToSmiles(frg)))
            if n_atoms == n_atoms_core:
                core_with_attachment = frg
            else:
                incomplete_frgs_core.append(frg)
                if (only_smallest_n_fragments is not None) and \
                   len(incomplete_frgs_core) >= only_smallest_n_fragments:
                    break
    # If fewer than only_smallest_n_fragments found, pad the list with None:
    else:
        if only_smallest_n_fragments is not None:
            incomplete_frgs_core += [None]*(only_smallest_n_fragments - len(incomplete_frgs_core))

    # If the core is not a BRICS fragment, optionally fragment further, back to the core:
    if not BRICS_only and (core_with_attachment is None):
        core_with_attachment = Chem.ReplaceSidechains(incomplete_frgs_core[0], core)

    return core_with_attachment, incomplete_frgs_core


def get_adjacent_brics_frgs(cmpd, substruct, BRICS_only=False):
    """
    Get BRICS fragments next to a particular substructure.
    """

    fn_grps_with_attachment = []
    fn_grps_capped = []
    fused_ring = []
    extended_core_frgs_smi = []

    # Convert substruct to mol object:
    if isinstance(substruct, str):
        substruct = Chem.MolFromSmiles(substruct)
   
    # Get molecule intermediates containing the substructure after BRICSDecompose:

    core_with_attachment, extended_core_frgs = \
    get_common_core_intermediates(cmpd, core=substruct, BRICS_only=BRICS_only)

    core_with_attachment = Chem.MolToSmiles(core_with_attachment)
    n_core_attachment_points = core_with_attachment.count('*')

    # Loop through intermediates:
    prev_n_frgs_attachment_points = 0
    for frg in extended_core_frgs:
        if frg is None:
            break

        # Cap open attachment points, remove substructure and count resulting 
        # attachment points:
        # Could be improved as involves lots of SMILES <-> Mol conversions:
        frg = Chem.MolFromSmiles(cap_attachment_points(Chem.MolToSmiles(frg)))
        adj_fn_grps = Chem.ReplaceCore(frg, 
                                       coreQuery=substruct)
        adj_fn_grps = Chem.MolToSmiles(adj_fn_grps) 
        n_frgs_attachment_points = adj_fn_grps.count('*')

        # Check how many of the core attachment points are accounted for, 
        # if this is greater than previous intermediates, then a new 
        # fragment must have been added at a separate attachment point:
        if n_frgs_attachment_points > prev_n_frgs_attachment_points:
            prev_n_frgs_attachment_points = n_frgs_attachment_points

            # Process fragments:
            fn_grps_with_attachment = []
            fn_grps_capped = []
            fused_ring = []
            extended_core_frgs_smi = []

            for adj_fn_grp in adj_fn_grps.split('.'):

                adj_fn_grp_capped = None
                canon_frg = None

                # Check if this fragment is attached by multiple points 
                # to the substructure (i.e. a fused ring):
                if adj_fn_grp.count('*') > 1:
                    fused_ring.append(True)
                else:
                    fused_ring.append(False)

                # If fused ring (and aromatic) RDKit will not read
                # it, but otherwise canonicalise:            
                if not fused_ring[-1]:
                    adj_fn_grp = canonicalise_frgms(adj_fn_grp, 
                                                    rm_BRICS_rules=True)
                    adj_fn_grp_capped = cap_attachment_points(adj_fn_grp)

                    for ori_frg in extended_core_frgs:
                        if ori_frg is None:
                            break
                        if ori_frg.HasSubstructMatch(Chem.MolFromSmiles(adj_fn_grp_capped)):
                            canon_frg = canonicalise_frgms(Chem.MolToSmiles(ori_frg), 
                                                           rm_BRICS_rules=True)
                            break

                fn_grps_with_attachment.append(adj_fn_grp)
                fn_grps_capped.append(adj_fn_grp_capped)
                extended_core_frgs_smi.append(canon_frg)

    return zip(fn_grps_with_attachment, fn_grps_capped, extended_core_frgs_smi, fused_ring)
