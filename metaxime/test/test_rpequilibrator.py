import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpEquilibrator

class TestRPEquilibrator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.rpeq = rpEquilibrator(path=os.path.join('data', 'rpequilibrator', 'rpsbml.xml'))
        cls.data = json.load(open(os.path.join('data', 'rpequilibrator', 'data.json'), 'r'))

    def test_makeSpeciesStr(self):
        self.assertEqual(self.rpeq._makeSpeciesStr(self.rpeq.model.getSpecies('MNXM89557__64__MNXC3')), 'CHEBI:5431')
        self.assertEqual(self.rpeq._makeSpeciesStr(self.rpeq.model.getSpecies('MNXM89557__64__MNXC3'), 'name'), 'L-glutamate')
        self.assertEqual(self.rpeq._makeSpeciesStr(self.rpeq.model.getSpecies('MNXM89557__64__MNXC3'), 'id'), 'MNXM89557')

    def test_makeReactionStr(self):
        self.assertEqual(self.rpeq._makeReactionStr(self.rpeq.model.getReaction('RP1')), '1.0 CHEBI:6280 + 1.0 CHEBI:5584 <=> 1.0 CHEBI:8650 + 1.0 CHEBI:3283 ')
        self.assertEqual(self.rpeq._makeReactionStr(self.rpeq.model.getReaction('RP1'), 'name'), '1.0 L-ornithine + 1.0 H(+) <=> 1.0 putrescine + 1.0 CO2 ')
        self.assertEqual(self.rpeq._makeReactionStr(self.rpeq.model.getReaction('RP1'), 'id'), '1.0 CMPD_0000000004 + 1.0 MNXM1 <=> 1.0 TARGET_0000000001 + 1.0 MNXM13 ')

    def test_speciesCmpQuery(self):
        res = self.rpeq._speciesCmpQuery(self.rpeq.model.getSpecies('MNXM89557__64__MNXC3'))
        self.assertEqual(res[0], self.data['speciescmpquery']['mu'])
        self.assertCountEqual(res[1], self.data['speciescmpquery']['sigma'])

    def test_reactionStrQuery(self):
        res = self.rpeq._reactionStrQuery(self.rpeq.model.getReaction('RP1'))
        self.assertTrue(res[0])
        self.assertAlmostEqual(res[1][0], -11.75708338789077)
        self.assertAlmostEqual(res[1][1], 1.4719330747677193)
        self.assertAlmostEqual(res[2][0], -67.87011607471129)
        self.assertAlmostEqual(res[2][1], 5.470350588406473)
        self.assertAlmostEqual(res[3][0], -26.579654512100888)
        self.assertAlmostEqual(res[3][1], 5.470350588406473)
        self.assertAlmostEqual(res[4][0], -43.69449204682192)
        self.assertAlmostEqual(res[4][1], 5.470350588406473)

    def test_pathway(self):
        res = self.rpeq.pathway(write_results=False)
        self.assertAlmostEqual(res[0][0], -53.0495385947695)
        self.assertAlmostEqual(res[0][1], 17.959120242335167)
        self.assertAlmostEqual(res[1][0], -87.27921366421157)
        self.assertAlmostEqual(res[1][1], 25.95074441357073)
        self.assertAlmostEqual(res[2][0], -53.0495385947695)
        self.assertAlmostEqual(res[2][1], 17.959120242335167)

    #def test_MDF(self) #<-- problem when the compounds are not found in the cc cache

    def test_toNetworkSBtab(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            hashlib.md5(open('/home/mdulac/Downloads/test.sbtab', 'rb').read()).hexdigest()
            self.assertTrue(self.rpeq.toNetworkSBtab(os.path.join(tmp_output_folder, 'test.sbtab')))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.sbtab'), 'rb').read()).hexdigest(), 'a0819ea10b15b670d8d2e0743674b0aa')

if __name__ == '__main__':
    unittest.main()
