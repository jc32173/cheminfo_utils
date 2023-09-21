# cheminfo\_utils

Set of useful chemoinformatics functions including curating/filtering datasets, manipulating/modifying molecules and calculating descriptors.  Most of the functions rely on RDKit, but some also use OpenEye which requires a license.

## Dataset curation and analysis

- `apply_lilly_rules.py`: Wrapper to call `Lilly_Medchem_Rules.rb` (https://github.com/IanAWatson/Lilly-Medchem-Rules) from python to apply Lilly rules to a list of SMILES and return the results in a DataFrame.

- `clustering.py`: Methods to cluster molecules by structures, including:
    - Butina clustering.

- `dataset_curation.py`: Set of functions for analysing a dataset of molecules, including calculating common representations and descriptors and identifying possible duplicates (including tautomers and stereoisomers).

## Processing molecules for ML

- `calc_desc.py`: Calculate molecular descriptors and fingerprints.

- `rdkit_extra_desc.py`: Extra molecular descriptors calculated using RDKit, including function to get number and size of fused ring systems (GetFusedRings).

- `smi_funcs.py`: Code for manipulating SMILES.
    - Canonicalising SMILES
    - Canonicalising tautomers
    - pH conversion using OpenEye or OpenBabel

## Modifying molecules

- `brics_fns.py`: Useful functions for fragmenting molecules using BRICS (see: https://chemistry-europe.onlinelibrary.wiley.com/doi/abs/10.1002/cmdc.200800178) and processing the resulting fragments.

- `join_mols.py`: Join two SMILES with attachment points with matching labels.

## Molecular structure

- `fn_grps.py`: Various functions for identifying functional groups next to a given substructure within molecules.  Includes Ertl's definition (https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0225-z).

- `struct_fns.py`: Code for calculating structural features of molecules, including:
    - The minimum through bond distance between atoms or substructures within a molecule.



