import unittest
from rdkit import Chem
import numpy as np
import sys
sys.path.append('/users/xpb20111/programs')
from rdkit_funcs import GetPossibleTautomers, canonicalise_tautomer


class Test_GetPossibleTautomers(unittest.TestCase):

    def setUp(self):
        """
        """

        self.smiles = 'Cc1cc2c(C3=CCOCC3)cnc(NC3CCN(CCC(N)=O)CC3)c2[nH]c1=O'
        self.mol = Chem.MolFromSmiles(self.smiles)
        self.canon_tauto = canonicalise_tautomer(self.smiles)

    def test_tautomer_oreder_1(self):
        """
        Check that tautomers are ordered from best to worst.
        """

        tautomers = GetPossibleTautomers(self.mol, 
                                         sort=True, 
                                         return_scores=False)
        self.assertEqual(self.canon_tauto, tautomers[0])
        
        #self.assertListEqual([-1, 4, 16, 4], n_bonds)

    def test_tautomer_order_2(self):
        """
        Check tautomer order.
        """

        tautomers = GetPossibleTautomers(self.mol, 
                                         sort=True, 
                                         return_scores=True)
        
        # NEED TO CHECK HOW TIES ARE RESOLVED
        tauto_sort_order = np.argsort([tauto[1] for tauto in tautomers])[::-1]
        print(list(zip(range(len(tautomers)), tautomers)))
        self.assertListEqual(list(tauto_sort_order), 
                             list(range(len(tautomers))))

    def test_tautomer_order_3(self):
        """
        Check whether tautomers are printed in order by default.
        """

        tautomers = GetPossibleTautomers(self.mol, 
                                         sort=False, 
                                         return_scores=True)
#        sorted_tautomers = GetPossibleTautomers(smiles, 
#                                                sort=True, 
#                                                return_scores=True)
        tauto_sort_order = np.argsort([tauto[1] for tauto in tautomers])
        #print(tauto_sort_order)
#        self.assertListEqual(self.canon_tauto, tautomers[0])
