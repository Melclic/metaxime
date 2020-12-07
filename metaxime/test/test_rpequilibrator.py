import unittest
import tarfile
import glob
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpEquilibrator
from metaxime import rpSBML

class TestRPEquilibrator(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.rpeq = rpEquilibrator(path=os.path.join('data', 'rpequilibrator', 'rpsbml.xml'))

    def test_runCollection(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpEquilibrator.runCollection(os.path.join('data', 'rpequilibrator', 'test.rpcol'),
                                         os.path.join(tmp_output_folder, 'test.rpcol'))
            tar = tarfile.open(os.path.join(tmp_output_folder, 'test.rpcol'), mode='r')
            os.mkdir(os.path.join(tmp_output_folder, 'results'))
            tar.extractall(os.path.join(tmp_output_folder, 'results'))
            self.assertTrue(len(glob.glob(os.path.join(tmp_output_folder, 'results', 'rpsbml_collection', 'models', '*')))==1)
            rpsbml = rpSBML(path=glob.glob(os.path.join(tmp_output_folder, 'results', 'rpsbml_collection', 'models', '*'))[0])
            asdict = rpsbml.asDict()
            self.assertAlmostEqual(asdict['pathway']['brsynth']['dfG_prime_o']['value'], 1784.7384959433493)

    def test_makeSpeciesStr(self):
        self.assertEqual(self.rpeq._makeSpeciesStr(self.rpeq.model.getSpecies('MNXM100__64__MNXC3')), 'CHEBI:5332')
        self.assertEqual(self.rpeq._makeSpeciesStr(self.rpeq.model.getSpecies('MNXM100__64__MNXC3'), 'name'), 'MNXM100')
        self.assertEqual(self.rpeq._makeSpeciesStr(self.rpeq.model.getSpecies('MNXM100__64__MNXC3'), 'id'), 'MNXM100__64__MNXC3')

    def test_makeReactionStr(self):
        self.assertEqual(self.rpeq._makeReactionStr(self.rpeq.model.getReaction('RP1')), '1.0 CHEBI:5332 <=> 1.0 XMGQYMWWDOXHJM-UHFFFAOYSA-N ')
        self.assertEqual(self.rpeq._makeReactionStr(self.rpeq.model.getReaction('RP1'), 'name'), '1.0 MNXM100 <=> 1.0 TARGET_0000000001 ')
        self.assertEqual(self.rpeq._makeReactionStr(self.rpeq.model.getReaction('RP1'), 'id'), '1.0 MNXM100__64__MNXC3 <=> 1.0 TARGET_0000000001__64__MNXC3 ')

    def test_speciesCmpQuery(self):
        mu, sigma = self.rpeq._speciesCmpQuery(self.rpeq.model.getSpecies('MNXM100__64__MNXC3'))
        self.assertAlmostEqual(mu, -1719.970141951422)

    def test_reactionStrQuery(self):
        res = self.rpeq._reactionStrQuery(self.rpeq.model.getReaction('RP1'))
        self.assertFalse(res[0])
        self.assertAlmostEqual(res[1][0], 720.3420273166071)
        self.assertAlmostEqual(res[1][1], 2.613560059758376)
        self.assertAlmostEqual(res[2][0], 1863.4700077198747)
        self.assertAlmostEqual(res[2][1], 6.475425663398934)
        self.assertAlmostEqual(res[3][0], 1784.7384959433493)
        self.assertAlmostEqual(res[3][1], 6.475425663398934)
        self.assertAlmostEqual(res[4][0], 1784.7384959433493)
        self.assertAlmostEqual(res[4][1], 6.475425663398934)

    def test_pathway(self):
        res = self.rpeq.pathway(write_results=False)
        self.assertAlmostEqual(res[0][0], 1784.7384959433493)
        self.assertAlmostEqual(res[0][1], 0.0)
        self.assertAlmostEqual(res[1][0], 1784.7384959433493)
        self.assertAlmostEqual(res[1][1], 0.0)
        self.assertAlmostEqual(res[2][0], 1784.7384959433493)
        self.assertAlmostEqual(res[2][1], 0.0)

    def test_MDF(self): 
        #TODO: find a file that does not fail
        rpeq = rpEquilibrator(path=os.path.join('data', 'rpequilibrator', 'rpsbml.xml'))
        res = rpeq.pathway(write_results=True)
        res = rpeq.MDF()
        self.assertAlmostEqual(res, 0.0)

    def test_toNetworkSBtab(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpeq = rpEquilibrator(path=os.path.join('data', 'rpequilibrator', 'rpsbml.xml'))
            res = rpeq.pathway(write_results=True)
            self.assertTrue(rpeq.toNetworkSBtab(os.path.join(tmp_output_folder, 'test.sbtab')))
            #TODO: read and and extract some values

if __name__ == '__main__':
    unittest.main()
