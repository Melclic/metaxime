import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpMerge

class TestRPMerge(unittest.TestCase):

    '''
    @classmethod
    def setUpClass(cls):
        #load a rpSBML file
        cls.data = json.load(open(os.path.join('data', 'rpmerge', 'data.json'), 'r'))
    '''

    def test_mergeSBMLFiles(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpMerge(os.path.join('data', 'rpmerge', 'rpsbml.xml'), os.path.join('data', 'rpmerge', 'gem.xml'), os.path.join(tmp_output_folder, 'test.sbml'))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.sbml'), 'rb').read()).hexdigest(), '66cb235c127e2bb07c6c13bea7bc6df2')

    #def _findUniqueRowColumn

    #def _checkSingleParent

    #def _compareReactions

    #def _containedReaction

    #def _compareReaction

    #def _compareSpecies

    #def _compareEC

    def test_mergeModels(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            gem = rpMerge('gem', path=os.path.join('data', 'rpmerge', 'gem.xml'))
            gem.mergeModels(os.path.join('data', 'rpmerge', 'rpsbml.xml'))
            gem.writeSBML(os.path.join(tmp_output_folder, 'test.sbml'))
            print(os.listdir(tmp_output_folder))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.sbml'), 'rb').read()).hexdigest(), '66cb235c127e2bb07c6c13bea7bc6df2')

if __name__ == '__main__':
    unittest.main()
