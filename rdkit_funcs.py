"""
Useful functions to read/manipulate molecules using RDKit mol objects and functions.
"""


import numpy as np
import re
import sys

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from openeye import oechem
from openeye import oequacpac

sys.path.insert(0, '/users/xpb20111/programs/deepchem_dev_nested_CV/cheminfo_utils')

from smi_funcs import adjust_for_ph

# Check OpenEye license:
if not oechem.OEChemIsLicensed():
    raise RuntimeError("No license found for OpenEye OEChem")
if not oequacpac.OEQuacPacIsLicensed():
    raise RuntimeError("No license found for OpenEye QuacPac")


#================#
# Load molecules #
#================#


def GetMol(smi):
    """
    Convert SMILES to RDKit mol object
    """
    mol = Chem.MolFromSmiles(smi)
    # Check smi can be read by rdkit:
    if mol is None:
        print('Error with smi: '+smi+', cannot convert to molecule using RDKit')
    return mol


# From: /users/xpb20111/Prosperity_partnership/Notebooks/202305_New_core/202305_New_core.ipynb
def get_mol_if_smiles(mol):
    """
    If function is given string as input, assume
    this is a SMILES and return an RDKit mol object.
    """
    
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            raise ValueError('Cannot convert {} to RDKit mol'.format(mol))
    return mol


# Check validaity:
def inchi_check(inchi):
    """
    Check InChI validity (used in Sorkun et al 2019).
    """
    mol_from_inchi = Chem.MolFromInchi(inchi)
    if mol_from_inchi is None:
        raise ValueError('Cannot convert InChI into RDKit Mol object.')
    if inchi != Chem.MolToInchi(mol_from_inchi):
        return False
    else:
        return True
    
    

#===================#
# Processing SMILES #
#===================#


# Functions for processing SMILES and curating datasets.

# Canonicalise SMILES:
def canonicalise_smiles(smi):
    """
    Convert SMILES to canonical RDKit form.
    """
    canon_smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi),
                                 canonical=True,
                                 isomericSmiles=True)
    return canon_smi


# Canonicalise tautomer:
def canonicalise_tautomer(smi, method='rdkit'):
    """
    Convert a SMILES to a canonical tautomer, 
    using RDKit, OpenBabel or via InChIs.
    """

    if method == 'rdkit':
        enumerator = rdMolStandardize.TautomerEnumerator()
        # enumerator.SetReassignStereo = True
        # enumerator.SetRemoveBondStereo = False
        # enumerator.SetRemoveSp3Stereo = False
        mol = enumerator.Canonicalize(Chem.MolFromSmiles(smi))
        canon_smi = Chem.MolToSmiles(mol, canonical=True,
                                     isomericSmiles=True)

        # If no change in SMILES, revert to original SMILES to retain
        # stereochemistry which may have been lost in 
        # TautomerEnumerator():
        canon_smi_2D = Chem.MolToSmiles(mol, canonical=True,
                                        isomericSmiles=False)
        smi_2D = Chem.MolToSmiles(Chem.MolFromSmiles(smi),
                                  canonical=True, isomericSmiles=False)
        if canon_smi != smi and canon_smi_2D == smi_2D:
            canon_smi = smi

    elif method == 'inchi':
        # Convert to and from InChI to standarise tautomer:
        mol_smi = Chem.MolFromSmiles(smi)
        inchi = Chem.MolToInchi(mol_smi)
        mol_inchi = Chem.MolFromInchi(inchi)
        canon_smi = Chem.MolToSmiles(mol_inchi, canonical=True,
                                     isomericSmiles=True)

    elif method == 'obabel':
        # Could use otautomer to get obabel canonical tautomers
        print('Warning: Method not yet implemented', file=sys.stderr)
        canon_smi = smi

    elif method == 'openeye':
        # Could use openeye to get a canonical tautomer
        print('Warning: Method not yet implemented', file=sys.stderr)
        canon_smi = smi

    return canon_smi


def GetPossibleTautomers(mol):
    """
    Generate list of all possible tautomers.
    """
    enumerator = rdMolStandardize.TautomerEnumerator()
    return [Chem.MolToSmiles(m) for m in enumerator.Enumerate(mol)]


# From: Troubleshoot_issues_with_acpype:
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
def GetPossibleStereoisomers(mol):
    """
    Generate list of all possible stereoisomers by enumerating over
    undefined stereochemistry.
    """
    opts = StereoEnumerationOptions(unique=True, onlyUnassigned=True)
    return [Chem.MolToSmiles(m) for m in EnumerateStereoisomers(mol,
                                                                options=opts)]


def GetPossibleTautomersStereoisomers(mol):
    smi_ls = []
    enumerator = rdMolStandardize.TautomerEnumerator()
    for tauto_mol in enumerator.Enumerate(mol):
        smi_ls += GetPossibleStereoisomers(tauto_mol)
    return smi_ls


def GetEnantiomer(smi):
    """
    Convert SMILES to enantiomer by reversing @ and @@ 
    labels for chiral sites.
    """
    if '@' in smi:
        smi = re.sub('@@', '__', smi)
        smi = re.sub('@', '@@', smi)
        smi = re.sub('__', '@', smi)
        # Check SMILES can be read:
        mol = GetMol(smi)
        if mol is None:
            return mol
        else:
            # Canonicalise:
            return Chem.MolToSmiles(mol, 
                                    canonical=True)
    else:
        return None


# From 202210_BRICS_hierachical_decomposition:
def get_frgm_fp(frgms, bit_smis):
    """
    Convert a set of fragments, such as from a BRICS decomposition
    of a molecule, to a fingerprint indicating which fragments from 
    a larger set are present.
    """
    fp = np.zeros(len(bit_smis), dtype=bool)
    for frg in frgms:
        fp[bit_smis.index(frg)] = 1
    return fp.astype(int).tolist()


# From:
def is_ionisable(smi, method='OpenEye'):
    """
    Test if SMILES is ionisable using OpenEye or OpenBabel.
    """

    if method == 'OpenEye':
        options = oequacpac.OEFormalChargeOptions()
        options.SetMaxCount(2)
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        if len(list(oequacpac.OEEnumerateFormalCharges(mol, options))) > 1:
            return True
        else:
            return False
    elif method == 'OpenBabel':
        if adjust_for_ph(smi, ph=1, phmodel='OpenBabel') != adjust_for_ph(smi, ph=14, phmodel='OpenBabel'):
            return True
        else:
            return False


# Function taken from /users/xpb20111/Solubility/AqSolDB/Dataset/Curate_dataset.ipynb
def ionised_smiles(smi, method='OpenEye'):
    smi_ls = []
    if method == 'OpenEye':
        options = oequacpac.OEFormalChargeOptions()
        options.SetMaxCount(20)
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        for mol in oequacpac.OEEnumerateFormalCharges(mol, options):
            smi_ls.append(oechem.OEMolToSmiles(mol))
    elif method == 'OpenBabel':
        pass
    return smi_ls


# From: /users/xpb20111/Prosperity_partnership/Notebooks/202305_New_core/202305_New_core.ipynb
def SmilesToSmartsScaffold(mol):
    """
    Convert SMILES to SMARTS scaffold, with all
    bonds replaced by any bond type.  Useful for
    substructure searching.
    """

    mol = Chem.MolFromSmiles(mol)
    sma = Chem.MolToSmarts(mol)
    sma = re.sub(':|-|=', '~', sma)
    return sma


#===============#
# Visualisation #
#===============#

# From: /users/xpb20111/Prosperity_partnership/Notebooks/202305_New_core/202305_New_core.ipynb
def save_mols_image(mols, legends=[], molsPerRow=4, out_filename=''):
    """
    Save grid of molecules to file as a png image.
    """

    im = \
    Chem.Draw.MolsToGridImage(mols=mols,
                              legends=legends,
                              molsPerRow=molsPerRow,
                              returnPNG=True)
    if out_filename != '':
        with open(out_filename, 'wb') as out:
            out.write(im.data)
        print('scp archie:'+os.path.abspath(out_file)+' .')
    try:
        display(im)
    except NameError:
        pass
    return im
