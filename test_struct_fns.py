import unittest
from struct_fns import calc_bond_distance_between_substructs


class Test_calc_bond_distance_between_substructs(unittest.TestCase):
    def test_1(self):
        n_bonds = calc_bond_distance_between_substructs('CCC(N)CCCCOCCCCC(=O)OCCCCc1ccccc1', 
                                                        'c1ccccc1', 
                                                        ['CNCN', 'COC', 'CN', 'COC(C)=O'])
        self.assertListEqual([-1, 4, 16, 4], n_bonds)
