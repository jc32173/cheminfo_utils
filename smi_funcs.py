# Writen October 2020

"""
Functions for modifying SMILES and calculating descriptors.
"""

__version__ = '2023.6.2'


from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
# Additional imports from rdkit for error logging:
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
# Need openbabel/openeye for some pH conversions:
from openbabel import openbabel as ob
import openeye as oe
from openeye import oechem
from io import StringIO
import sys
import os
import numpy as np
import pandas as pd
import re
import logging

from cheminfo_utils.rdkit_extra_desc import extra_rdkit_descs


# Set up logger for module:
logger = logging.getLogger(__name__)

# Silence RDKit warnings:
logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog('rdApp.error')


# Canonicalise SMILES:
def canonicalise_smiles(smi, method='RDKit'):
    """
    Convert SMILES to canonical RDKit form.
    """

    canon_smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi), 
                                 canonical=True,
                                 isomericSmiles=True)
    return canon_smi


# Canonicalise tautomer:
def canonicalise_tautomer(smi, method='RDKit'):
    """
    Convert a SMILES to a canonical tautomer, 
    using RDKit, OpenBabel or via InChIs.
    """

    # Access any rdkit errors from stderr:
    #Chem.WrapLogs()
    sio = sys.stderr = StringIO()

    if method == 'RDKit':
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

    elif method == 'InChI':
        # Convert to and from InChI to standarise tautomer:
        mol_smi = Chem.MolFromSmiles(smi)
        inchi = Chem.MolToInchi(mol_smi)
        mol_inchi = Chem.MolFromInchi(inchi)
        canon_smi = Chem.MolToSmiles(mol_inchi, canonical=True, 
                                     isomericSmiles=True)

    elif method == 'OpenBabel':
        # Could use otautomer to get obabel canonical tautomers
        print('Warning: Method not yet implemented', file=sys.stderr)
        canon_smi = smi

    elif method == 'OpenEye':
        # Could use openeye to get a canonical tautomer
        print('Warning: Method not yet implemented', file=sys.stderr)
        canon_smi = smi

    # Get any warnings in stderr:
    warnings = sio.getvalue()

    # Redirect errors back to stderr:
    sys.stderr = sys.__stderr__

    return canon_smi, warnings


# Canonicalise tautomer, but preserve given substructure tautomer form:
def canonicalise_tautomer_substruct(cmpd, substruct=None):
    """
    Select the canonical tautomer as the highest scoring tautomer which has a 
    given tautomeric form for a common substructure.
    
    >>> Chem.MolToSmiles(
    ... canonicalise_tautomer_substruct(cmpd='CC(=O)CCCCCC(=O)CC', 
    ...                                 substruct='C=C(O)C'))
    'CCC(O)=CCCCCC(C)=O'
    >>> Chem.MolToSmiles(
    ... canonicalise_tautomer_substruct(
    ...     cmpd='C2[C-]([O-])[N+](C=[N+](O)O)=C([O-])N2C[N+](C)(C)O',
    ...     substruct='C1NC(=O)NC(=O)1'))
    WARNING: No tautomers found with given substructure form, will return canonical tautomer.
    'C[N+](C)(O)CN1C[C-]([O-])[N+](C=[N+](O)O)=C1[O-]'
    """
    
    enumerator = rdMolStandardize.TautomerEnumerator()
    
    if isinstance(cmpd, str):
        cmpd = Chem.MolFromSmiles(cmpd)

    if substruct is None:
        return enumerator.Canonicalize(cmpd)
        
    if isinstance(substruct, str):
        substruct = Chem.MolFromSmiles(substruct)

    tauto_ls = []
    for i, tauto in enumerate(enumerator.Enumerate(cmpd)):
        tauto_ls.append([enumerator.ScoreTautomer(tauto), tauto])

    # Sort based on tautomer scores:
    tauto_ls = sorted(tauto_ls, key=lambda x: x[0])[::-1]

    for _, tauto in tauto_ls:
        if tauto.HasSubstructMatch(substruct):
            return tauto

    print('WARNING: No tautomers found with given substructure form, will return canonical tautomer.')
    return enumerator.Canonicalize(cmpd)


#def select_tautomer(rank=0):
#    """
#    Select specific ranked tautomer or one tautomer at random.
#    """
#    for i, tauto in enumerate(enumerator.Enumerate(cmpd)):
#        tauto_ls.append([enumerator.ScoreTautomer(tauto), tauto])


# Correct SMILES using regex:
def correct_smiles(smi, smi_transforms, check_rdkit=True):
    """
    Modify SMILES using regex.
    """

    # Access any rdkit errors from stderr:
    #Chem.WrapLogs()
    sio = sys.stderr = StringIO()

    for match_pat, repl_pat in smi_transforms:
        corr_smi = re.sub(match_pat, repl_pat, smi)

    # Check that new SMILES can be loaded into RDKit:
    if check_rdkit == True and Chem.MolFromSmiles(corr_smi) is None:
        warnings = sio.getvalue() + \
        'WARNING (correct_smiles): '
        'Error during tautomer correction, '
        'will use original SMILES.\n'
        corr_smi = smi
    else:
        warnings = sio.getvalue()

    # Redirect errors back to stderr:
    sys.stderr = sys.__stderr__

    return corr_smi, warnings


# Protonate/deprotonate to get SMILES at a given pH:
def adjust_for_ph(smi, 
                  ph=7.4, 
                  phmodel='OpenEye', 
                  phmodel_dir=None, 
                  #verbose=0 # False
                 ):
    """
    Protonate/deprotonate SMILES according
    to a given pH using OpenBabel or OpenEye.
    """

    if phmodel == 'OpenEye' and ph != 7.4:
        raise ValueError('Cannot use OpenEye pH conversion for pH != 7.4')

    #if verbose:
    #    print('Adjusting SMILES to pH: {}, using {}'.format(ph, phmodel))
    logger.info('Adjusting SMILES to pH: {}, using {}'.format(ph, phmodel))

    if phmodel == 'OpenBabel':

        # Access any rdkit errors from stderr:
        #Chem.WrapLogs()
        sio = sys.stderr = StringIO()

        # Set BABEL_DATADIR environment variable if using a modified 
        # phmodel.txt file containing pH transformations not in the
        # current directory
        if phmodel_dir is not None:
            os.environ['BABEL_DATADIR'] = phmodel_dir

        # Use obabel from python to do pH correction:
        ob_conv = ob.OBConversion()
        ob_mol = ob.OBMol()
        ob_conv.SetInAndOutFormats("smi", "smi")
        ob_conv.ReadString(ob_mol, smi)
        ob_mol.AddHydrogens(False,  # only add polar H (i.e. not to C atoms)
                            True,   # correct for pH
                            ph)     # pH value
        ph_smi = ob_conv.WriteString(ob_mol,
                                     True)  # trimWhitespace

        # Check that pH adjusted SMILES can be read by RDKit, 
        # if not return original SMILES and a warning:
        if Chem.MolFromSmiles(ph_smi) is None:
            warnings = sio.getvalue() + \
                'WARNING (adjust_for_ph): ' \
                'Error during pH correction with obabel, ' \
                'will use original SMILES.\n'
            ph_smi = smi
        else:
            warnings = sio.getvalue()

        # Redirect errors back to stderr:
        sys.stderr = sys.__stderr__

    elif phmodel == 'OpenEye':

        # Save any OE errors:
        warnos = oechem.oeosstream()
        oechem.OEThrow.SetOutputStream(warnos)

        # Need to add OE error handling, OEThrow?
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, smi)
        oe.OESetNeutralpHModel(mol)
        ph_smi = oechem.OEMolToSmiles(mol)

        # Get warnings from output stream, gives a
        # byte-string, so also have to decode to Unicode:
        warnings = str(warnos.str(), 'UTF-8')

        # Redirect OE errors back to stderr:
        oechem.OEThrow.SetOutputStream(oe.oeerr)

    return ph_smi, warnings


# Process SMILES before calculating descriptors
# (canonicalise tautomer, correct smiles, adjust for ph)
def process_smiles(smi, tauto=False, ph=None, phmodel=None, canon_smiles=False):
    """
    Process SMILES by converting to canonical tautomer, specific pH and 
    canonical SMILES.
    """
    warnings = ''

    if tauto:
        # Convert to canonical tautomer using rdkit:
        smi, warning = canonicalise_tautomer(smi, method='RDKit')
        if warning != '':
            warnings += warning

    # Correct known errors in SMILES from canonicalizeTautomers in rdkit:
    smi, warning = correct_smiles(smi, [[r"\[PH\]\(=O\)\(=O\)O",
                                         r"P(=O)(O)O"]])
    if warning != '':
        warnings += warning

    # Change tetrazole tautomer by moving H to first N (i.e. c1[nH]nnn1)
    # so that it can be correctly deprotonated in obabel if pH < pKa:
    smi, warning = correct_smiles(smi, [[r"c(\d+)nn\[nH\]n(\d+)",
                                         r"c\1[nH]nnn\2"]])
    if warning != '':
        warnings += warning

    if ph is not None:
        # Convert SMILES to given pH:
        smi, warning = adjust_for_ph(smi, ph=ph, phmodel=phmodel)
        if warning != '':
            warnings += warning

    # Make SMILES canonical:
    if canon_smiles:
        smi = canonicalise_smiles(smi)

    return smi, warnings


# Calculate rdkit descriptors:
def calc_rdkit_descs(cmpd_ls, 
                     id_ls=[], 
                     #input_type='SMILES', 
                     desc_ls=[], 
                     extra_desc=[]):
    """
    Calculate RDKit descriptors from SMILES.
    """

    if isinstance(cmpd_ls, str):
        cmpd_ls = [cmpd_ls]

    # Set up index:
    if len(id_ls) == 0:
        if isinstance(cmpd_ls[0], str):
            id_ls = cmpd_ls
        else:
            id_ls = range(len(cmpd_ls))

    # Combine index and compounds:
    cmpd_id_ls = zip(id_ls, cmpd_ls)

    ## Access any rdkit errors from stderr:
    #Chem.WrapLogs()
    #sio = sys.stderr = StringIO()

    # If no descriptors given, calculate all available RKDit 
    # descriptors:
    if len(desc_ls) == 0:
        desc_ls = [x[0] for x in Chem.Descriptors._descList]

    # Calculate descriptors in desc_ls:
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(desc_ls)
    df_descs = pd.DataFrame(data=np.zeros((len(cmpd_ls), len(desc_ls)+len(extra_desc))),
                            index=id_ls,
                            columns=desc_ls+extra_desc)
    df_descs.loc[:,:] = np.nan

    for cmpd_i, cmpd in cmpd_id_ls:

        # If string, assume this is a SMILES:
        if isinstance(cmpd, str):
            mol = Chem.MolFromSmiles(smi)
            # Check SMILES has been read by rdkit:
            if mol is None:
                raise ValueError('Cannot read SMILES: {} using RDKit.'.format(smi))

        # Otherwise assume this is an RDKit Mol object:
        else:
            mol = cmpd
        if mol is None:
            raise ValueError('Mol is None.')

        df_descs.loc[cmpd_i, desc_ls] = np.array(calc.CalcDescriptors(mol))

        # Calculate Ipc descriptor separately to ensure average value is
        # calculated, this is not default in RDKit-2020.09.1, but otherwise 
        # leads to very large values which can cause problems:
        if 'Ipc' in desc_ls:
            df_descs.loc[cmpd_i, 'Ipc'] = Chem.GraphDescriptors.Ipc(mol, 
                                                                    avg=True)
        for desc_name, desc_fn in extra_rdkit_descs:
            if desc_name in extra_desc:
                df_descs.loc[cmpd_i, desc_name] = desc_fn(mol)

    warnings = '' #sio.getvalue()

    ## Redirect errors back to stderr:
    #sys.stderr = sys.__stderr__

    return df_descs, warnings
