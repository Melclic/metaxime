import unittest
import sys
import os
import json
import hashlib
import tempfile
import libsbml

sys.path.insert(0, '../..')

from metaxime import rpCache
from metaxime import rpReader

class TestRPreader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data = json.load(open(os.path.join('data', 'rpreader', 'data.json'), 'r'))
        cls.rpcache = rpCache()
        cls.rpcache.populateCache()
        cls.rpreader = rpReader(self.rpcache)

    def test_rp2ToCollection(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpReader.rp2ToCollection(os.path.join('data', 'rpreader', 'rp_pathways.csv'),
                                     os.path.join('data', 'rpreader', 'rp2paths_compounds.csv'),
                                     os.path.join('data', 'rpreader', 'rp2paths_pathways.csv'),
                                     os.path.join(tmp_output_folder, 'test.rpcol'))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.rpcol'), 'rb').read()).hexdigest(), '4efe1038eb2c0ccc12e04f44dfc8e1a0')

    def test_rp2ToSBML(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            os.mkdir(os.path.join(tmp_output_folder, 'rpsbml_collection'))
            os.mkdir(os.path.join(tmp_output_folder, 'rpsbml_collection', 'models'))
            status = self.rpreader.rp2ToSBML(os.path.join('data', 'rpreader', 'rp_pathways.csv'),
                                             os.path.join('data', 'rpreader', 'rp2paths_compounds.csv'),
                                             os.path.join('data', 'rpreader', 'rp2paths_pathways.csv'),
                                             os.path.join(tmp_output_folder, 'rpsbml_collection', 'models'))
            with tarfile.open(os.path.join(tmp_output_folder, 'test.rpcol'), "w:xz") as tar:
                tar.add(os.path.join(tmp_output_folder, 'rpsbml_collection'), arcname='rpsbml_collection')
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.rpcol'), 'rb').read()).hexdigest(), '4efe1038eb2c0ccc12e04f44dfc8e1a0')

    def test_readRp2PathsCompounds(self):
        self.assertDictEqual(self.rpreader.readRp2PathsCompounds(os.path.join('data', 'rpreader', 'rp2paths_compounds.csv')), data['readrp2pathscompounds'])

    def test_readRpPathways(self):
        self.assertDictEqual(self.rpreader.readRpPathways(os.path.join('data', 'rpreader', 'rp_pathways.csv')), data['readrppathways'])

    def test_readRp2PathsPathways(self):
        self.assertDictEqual(self.rpreader.readRp2PathsPathways(os.path.join('data', 'rpreader', 'rp2paths_pathways.csv')), data['readrp2pathspathways'])

    #def test_gatherSpeciesInfo(self):

    #def test_completeReac(self):

    #def test_addCofactorSpecies(self):

    #def test_updateStoichio(self):

    #def test_addCofactorsStep(self):

    #def test_addCofactors(self):

    #def test_parseTSV(self):

    #def test_TSVtoSBML(self):
