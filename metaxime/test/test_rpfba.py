import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys
import tarfile
import glob

sys.path.insert(0, '../..')

from metaxime import rpFBA
from metaxime import rpSBML

class TestRPFBA(unittest.TestCase):

    """
    @classmethod
    def setUpClass(cls):
        #load a rpSBML file
        #self.data = json.load(open(os.path.join('data', 'rpfba', 'data.json'), 'r'))
        pass
    """

    def test_runCollection(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpFBA.runCollection(os.path.join('data', 'rpfba', 'test.rpcol'),
                                os.path.join('data', 'rpfba', 'gem.xml'),
                                os.path.join(tmp_output_folder, 'test.rpcol'),
                                num_workers=1)
            tar = tarfile.open(os.path.join(tmp_output_folder, 'test.rpcol'), mode='r')
            os.mkdir(os.path.join(tmp_output_folder, 'results'))
            tar.extractall(os.path.join(tmp_output_folder, 'results'))
            self.assertTrue(len(glob.glob(os.path.join(tmp_output_folder, 'results', 'rpsbml_collection', 'models', '*')))==1)
            rpsbml = rpSBML(path=glob.glob(os.path.join(tmp_output_folder, 'results', 'rpsbml_collection', 'models', '*'))[0])
            asdict = rpsbml.asDict()
            self.assertAlmostEqual(asdict['pathway']['brsynth']['fba_obj_biomass_restricted']['value'], 0.6577479108178588)
            self.assertAlmostEqual(asdict['pathway']['brsynth']['fba_obj_fraction']['value'], 0.9438866396863238)
            self.assertAlmostEqual(asdict['pathway']['brsynth']['fba_obj_biomass']['value'], 0.8769972144238116)

    #def _convertToCobra

    #def _writeAnalysisResults

    #def runMultiObjective

    def test_runFBA(self):
        rpfba = rpFBA(sbml_path=os.path.join('data', 'rpfba', 'gem.xml'))
        rpfba.mergeModels(os.path.join('data', 'rpfba', 'rpsbml.xml'))
        rpfba.runFBA('biomass', 1.0, write_results=True)
        asdict = rpfba.asDict()
        self.assertAlmostEqual(asdict['pathway']['brsynth']['fba_obj_biomass']['value'], 0.8769972144238116)

    def test_runParsimoniousFBA(self):
        rpfba = rpFBA(sbml_path=os.path.join('data', 'rpfba', 'gem.xml'))
        rpfba.mergeModels(os.path.join('data', 'rpfba', 'rpsbml.xml'))
        rpfba.runParsimoniousFBA('biomass', 1.0, write_results=True)
        asdict = rpfba.asDict()
        self.assertAlmostEqual(asdict['pathway']['brsynth']['fba_obj_biomass']['value'], 732.4143562781677)

    def test_runFractionReaction(self):
        #rpfba = rpFBA(rpsbml_path=os.path.join('data', 'rpfba', 'rpsbml.xml'), gem_sbml_path=os.path.join('data', 'rpfba', 'gem.xml'))
        rpfba = rpFBA(sbml_path=os.path.join('data', 'rpfba', 'gem.xml'))
        rpfba.mergeModels(os.path.join('data', 'rpfba', 'rpsbml.xml'))
        rpfba.runFractionReaction('biomass', 1.0, 'RP1_sink', 1.0, write_results=True)
        asdict = rpfba.asDict()
        self.assertAlmostEqual(asdict['pathway']['brsynth']['fba_obj_biomass_restricted']['value'], 0.6577479108178588)
        self.assertAlmostEqual(asdict['pathway']['brsynth']['fba_obj_fraction']['value'], 0.9438866396863238)
        self.assertAlmostEqual(asdict['pathway']['brsynth']['fba_obj_biomass']['value'], 0.8769972144238116)

if __name__ == '__main__':
    unittest.main()
