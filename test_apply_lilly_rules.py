import unittest
import pandas as pd
from pandas.testing import assert_frame_equal
from cheminfo_utils.apply_lilly_rules import apply_lilly_rules

class Test_apply_lilly_rules(unittest.TestCase):
    def test_1(self, 
               lilly_rules_script=\
        '/users/xpb20111/software/Lilly-Medchem-Rules/Lilly_Medchem_Rules.rb'):

        expected_results = \
        pd.DataFrame(data=[['CCCCCCC(=O)O', 'CCCCCCC(=O)O', True, 
                            'D(80) C6:no_rings', 'CCCCCCC(=O)O'], 
                           ['CCC', 'CCC', False, 
                            'TP1 not_enough_atoms', 'CCC'], 
                           ['CCCCC(=O)OCC', 'CCCCC(=O)OCC', True, 
                            'D(75) ester:no_rings:C4', 'CCCCC(=O)OCC'], 
                           ['c1ccccc1CC(=O)C', 'CC(=O)CC1=CC=CC=C1', True, 
                            None, 'CC(=O)CC1=CC=CC=C1']], 
                     columns=['SMILES', 'SMILES_Kekule', 'Lilly_rules_pass', 
                              'Lilly_rules_warning', 'Lilly_rules_SMILES'], 
                     index=['0', '1', '2', '3'])
        
        test_results = apply_lilly_rules(smiles=['CCCCCCC(=O)O', 
                                                 'CCC', 
                                                 'CCCCC(=O)OCC', 
                                                 'c1ccccc1CC(=O)C'], 
                                         lilly_rules_script=lilly_rules_script)

        assert_frame_equal(expected_results, test_results)
