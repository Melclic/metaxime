import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpFBA

class TestRPFBA(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        #load a rpSBML file
        #self.data = json.load(open(os.path.join('data', 'rpfba', 'data.json'), 'r'))
        pass

    #def _convertToCobra

    #def _writeAnalysisResults

    #def mergeModels

    def test_writeSBML(self):
        rpfba = rpFBA(rpsbml_path=os.path.join('data', 'rpfba', 'rpsbml.xml'), gem_sbml_path=os.path.join('data', 'rpfba', 'gem.xml'))
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpfba.writeSBML(os.path.join(tmp_output_folder, 'merged.xml'))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'merged.xml'), 'rb').read()).hexdigest(), '0e755a7ae4605279df728b5dab176181')

    #def runMultiObjective

    #def runFBA

    #def runParsimoniousFBA

    def test_runFractionReaction():
        rpfba = rpFBA(rpsbml_path=os.path.join('data', 'rpfba', 'rpsbml.xml'), gem_sbml_path=os.path.join('data', 'rpfba', 'gem.xml'))
        rpfba.runFractionReaction('biomass', 1.0, 'RP1_sink', 1.0)
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpfba.writeSBML(os.path.join(tmp_output_folder, 'merged.xml'))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'merged.xml'), 'rb').read()).hexdigest(), '0e755a7ae4605279df728b5dab176181')

if __name__ == '__main__':
    unittest.main()
