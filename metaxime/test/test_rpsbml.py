import unittest
import os
import json

from metaxime import rpSBML

class TestRPSBML(unittest.TestCase):
    def setUp(self):
        #load a rpSBML file
        self.rpsbml = rpSBML('test', path=os.path.join(os.path.dirname(__file__), 'data', 'rpsbml.xml'))
        self.data = json.load(os.path.join(os.path.dirname(__file__), 'data', 'data.json'))

    def test_computeMeanRulesScore(self):
        self.assertAlmostEqual(self.rpsbml._computeMeanRulesScore(), 0.5684564101634014, places=7, message='Equal to 0.5684564101634014')

    def test_dictRPpathway(self):
        self.assertDictEqual(self.rpsbml._dictRPpathway(), self.data['dictrppathway'])

    def test_nameToSbmlId(self):
        self.assertEqual(self.rpsbml._nameToSbmlId('test123-_!"£$%^&*(){}@~><>?'), 'test123___________________')

    def test_genMetaID(self):
        self.assertEqual(self.rpsbml._genMetaID('test123'), 'cc03e747a6afbbcbf8be7668acfebee5')

    def test_asDict(self):
        self.assertDictEqual(self.rpsbml.asDict(), self.data['asdict'])

    def test_readRPrules(self):
        self.assertDictEqual(self.rpsbml.readRPrules(), self.data['readrprules'])

    def test_getGroupsMembers(self):
        self.assertCountEqual(self.rpsbml.getGroupsMembers('rp_pathway'), ['RP1', 'RP2', 'RP3'])

    def test_readRPspecies(self):
        self.assertDictEqual(self.rpsbml.readRPspecies(), self.data['readrpspecies'])

    def test_readUniqueRPspecies(self):
        self.assertCountEqual(self.rpsbml.readUniqueRPspecies(), ['TARGET_0000000001__64__MNXC3', 'MNXM13__64__MNXC3', 'CMPD_0000000004__64__MNXC3', 'MNXM1__64__MNXC3', 'MNXM20__64__MNXC3', 'CMPD_0000000013__64__MNXC3', 'MNXM89557__64__MNXC3', 'MNXM5__64__MNXC3', 'MNXM7__64__MNXC3', 'MNXM9__64__MNXC3', 'MNXM6__64__MNXC3', 'MNXM3__64__MNXC3'])

if __name__ == '__main__':
    unittest.main()
