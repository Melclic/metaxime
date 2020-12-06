import unittest
import sys
import os
import json
import hashlib
import tempfile
import tarfile
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
        cls.rpreader = rpReader(rpcache=cls.rpcache)
        cls.maxDiff = None

    def test_rp2ToCollection(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpReader.rp2ToCollection(os.path.join('data', 'rpreader', 'rp_pathways.csv'),
                                     os.path.join('data', 'rpreader', 'rp2paths_compounds.csv'),
                                     os.path.join('data', 'rpreader', 'rp2paths_pathways.csv'),
                                     os.path.join(tmp_output_folder, 'test.rpcol'))
            tar = tarfile.open(os.path.join(tmp_output_folder, 'test.rpcol'), mode='r')
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.rpcol'), 'rb').read()).hexdigest(), 'b3ae498ea2b5dbcc2c0e35410b875865')

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
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.rpcol'), 'rb').read()).hexdigest(), '4d802a41c3bd6fc4bdd53890c9de5a7d')

    def test_readRp2PathsCompounds(self):
        self.assertDictEqual(self.rpreader.readRp2PathsCompounds(os.path.join('data', 'rpreader', 'rp2paths_compounds.csv')), self.data['readrp2pathscompounds'])

    '''
    def test_readRpPathways(self):
        self.assertDictEqual(self.rpreader.readRpPathways(os.path.join('data', 'rpreader', 'rp_pathways.csv')), self.data['readrppathways'])
        rp_transformations, rp_compounds = (self.rpreader.readRpPathways(os.path.join('data', 'rpreader', 'rp_pathways.csv'))
        self.assertDictEqual(rp_transformations, self.data['readrppathways']['rp_transformations'])
        self.assertDictEqual(rp_compounds, self.data['readrppathways']['rp_compounds'])
    '''

    def test_readRp2PathsPathways(self):
        self.assertDictEqual(self.rpreader.readRp2PathsPathways(os.path.join('data', 'rpreader', 'rp2paths_pathways.csv')), self.data['readrp2pathspathways'])

    #def test_gatherSpeciesInfo(self):

    #def test_completeReac(self):

    #def test_addCofactorSpecies(self):

    #def test_updateStoichio(self):

    #def test_addCofactorsStep(self):

    #def test_addCofactors(self):

    #def test_parseTSV(self):

    #def test_TSVtoSBML(self):

if __name__ == '__main__':
    unittest.main()
