import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpMerge
from metaxime import rpSBML

class TestRPMerge(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.gem = rpSBML(path=os.path.join('data', 'rpmerge', 'gem.xml'))
        self.rpsbml = rpSBML(path=os.path.join('data', 'rpmerge', 'rpsbml.xml'))
        self.rpmerge = rpMerge(True, path=os.path.join('data', 'rpmerge', 'gem.xml'))

    def test_mergeSBMLFiles(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpMerge.mergeSBMLFiles(os.path.join('data', 'rpmerge', 'rpsbml.xml'),
                                   os.path.join('data', 'rpmerge', 'gem.xml'),
                                   os.path.join(tmp_output_folder, 'merged.xml'))
            rpsbml = rpSBML(path=os.path.join(tmp_output_folder, 'merged.xml'))
            reac = rpsbml.model.getReaction('RP1')
            self.assertEqual(reac.getId(), 'RP1')
            self.assertTrue(len([i.species for i in reac.getListOfReactants()])==1)
            spe_id = [i.species for i in reac.getListOfReactants()][0]
            self.assertEqual(spe_id, 'M_grdp_c')
            self.assertTrue(len([i.species for i in reac.getListOfProducts()])==1)
            spe_id = [i.species for i in reac.getListOfProducts()][0]
            self.assertEqual(spe_id, 'TARGET_0000000001__64__MNXC3')

    #def _findUniqueRowColumn

    #def _checkSingleParent

    def test_compareReactions(self):
        species_source_target = self.rpmerge._compareSpecies({'MNXC3': 'MNXC3', 'MNXC2': 'MNXC2', 'MNXC19': 'MNXC19'}, self.rpsbml.model, self.gem.model)
        match = self.rpmerge._compareReactions(species_source_target, self.rpsbml.model, self.gem.model)
        self.assertDictEqual(match, {'RP1': {}, 'RP1_sink': {}})

    def test_compareReaction(self):
        species_source_target = self.rpmerge._compareSpecies({'MNXC3': 'MNXC3', 'MNXC2': 'MNXC2', 'MNXC19': 'MNXC19'}, self.rpsbml.model, self.gem.model)
        match, score = self.rpmerge._compareReaction(species_source_target, self.rpsbml.model.getReaction('RP1'), self.gem.model.getReaction('R_MGSA'), ec_match=True) 
        self.assertFalse(match)
        self.assertAlmostEqual(score, 0.1875)

    def test_compareSpecies(self):
        match = self.rpmerge._compareSpecies({'MNXC3': 'MNXC3', 'MNXC2': 'MNXC2', 'MNXC19': 'MNXC19'}, self.rpsbml.model, self.gem.model)
        self.assertDictEqual(match, {'MNXM100__64__MNXC3': {'M_grdp_c': 0.4}, 'TARGET_0000000001__64__MNXC3': {}})

    def test_compareSingleSpecies(self):
        e = self.rpsbml.model.getSpecies('MNXM100__64__MNXC3')
        f = self.gem.model.getSpecies('M_dhap_c')
        match, score = self.rpmerge._compareSingleSpecies(e,f)
        self.assertFalse(match)
        self.assertAlmostEqual(score, 0.0)
        match, score = self.rpmerge._compareSingleSpecies(e,e)
        self.assertTrue(match)
        self.assertAlmostEqual(score, 1.0)

    def test_compareEC(self):
        self.assertAlmostEqual(self.rpmerge._compareEC(self.rpsbml.model.getReaction('RP1'), self.gem.model.getReaction('R_MGSA')), 0.75)

    def test_mergeModels(self):
        gem = rpMerge(True, path=os.path.join('data', 'rpmerge', 'gem.xml'))
        #check to make sure that indeed the reaction does not exist
        self.assertIsNone(gem.model.getReaction('RP1'))
        gem.mergeModels(os.path.join('data', 'rpmerge', 'rpsbml.xml'))
        #Retreive a rpsbml reaction and make sure that the new species have been created and the existing ones are there
        reac = gem.model.getReaction('RP1')
        self.assertEqual(reac.getId(), 'RP1')
        self.assertTrue(len([i.species for i in reac.getListOfReactants()])==1)
        spe_id = [i.species for i in reac.getListOfReactants()][0]
        self.assertEqual(spe_id, 'M_grdp_c')
        self.assertTrue(len([i.species for i in reac.getListOfProducts()])==1)
        spe_id = [i.species for i in reac.getListOfProducts()][0]
        self.assertEqual(spe_id, 'TARGET_0000000001__64__MNXC3')

if __name__ == '__main__':
    unittest.main()
