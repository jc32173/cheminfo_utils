from rdkit import Chem
import pandas as pd
import os


def mols_to_grid(cmpds, 
                 legends=None, 
                 molsPerRow=5, 
                 filename=False, 
                 notebook=True):
    """
    Draw molecule structures in a grid using RDKit MolsToGridImage and 
    display the results on a notebook and/or save them to a png file.
    """

    if isinstance(cmpds, str):
        cmpds = [Chem.MolFromSmiles(cmpds)]
    elif isinstance(cmpds, pd.core.series.Series):
        cmpds = cmpds.to_list()
    elif isinstance(cmpds, pd.core.frame.DataFrame):
        cmpds = cmpds.squeeze().to_list()
        
    if isinstance(cmpds[0], str):
        cmpds = [Chem.MolFromSmiles(cmpd) for cmpd in cmpds]

    im = Chem.Draw.MolsToGridImage(mols=cmpds, 
                                   legends=legends, 
                                   molsPerRow=molsPerRow)
    if notebook:
        display(im)

    if filename:
        with open(filename, 'wb') as out:
            out.write(im.data)
        print("scp archie:"+os.getcwd()+"/"+filename+" .")

    return im


def save_mol_to_file():
    pass
