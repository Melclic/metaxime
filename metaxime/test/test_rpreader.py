import unittest
import sys
import os
import json
import tempfile
import glob
import tarfile
import libsbml

sys.path.insert(0, '../..')

from metaxime import rpCache
from metaxime import rpReader
from metaxime import rpSBML

class TestRPreader(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.data = json.load(open(os.path.join('data', 'rpreader', 'data.json'), 'r'))
        self.rpreader = rpReader()
        self.maxDiff = None

    def test_rp2ToCollection(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpReader.rp2ToCollection(os.path.join('data', 'rpreader', 'rp_pathways.csv'),
                                     os.path.join('data', 'rpreader', 'rp2paths_compounds.tsv'),
                                     os.path.join('data', 'rpreader', 'rp2paths_pathways.csv'),
                                     os.path.join(tmp_output_folder, 'test.rpcol'))
            tar = tarfile.open(os.path.join(tmp_output_folder, 'test.rpcol'), mode='r')
            os.mkdir(os.path.join(tmp_output_folder, 'results'))
            tar.extractall(os.path.join(tmp_output_folder, 'results'))
            #make sure that rpreader only creates a single model file
            self.assertTrue(len(glob.glob(os.path.join(tmp_output_folder, 'results', 'rpsbml_collection', 'models', '*')))==1)
            #compare the results to the data.json results
            rpsbml = rpSBML(path=glob.glob(os.path.join(tmp_output_folder, 'results', 'rpsbml_collection', 'models', '*'))[0])
            #WARNING: due to the nature of the json writting and reading, that are all strings, we perform the same operation here
            json.dump({'asDict': rpsbml.asDict()}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['asDict'], self.data['asDict'])

    def test_rp2ToSBML(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            status = self.rpreader.rp2ToSBML(os.path.join('data', 'rpreader', 'rp_pathways.csv'),
                                             os.path.join('data', 'rpreader', 'rp2paths_compounds.tsv'),
                                             os.path.join('data', 'rpreader', 'rp2paths_pathways.csv'),
                                             os.path.join(tmp_output_folder))
            #make sure that rpreader only creates a single model file
            self.assertTrue(len(glob.glob(os.path.join(tmp_output_folder, '*')))==1)
            #compare the results to the data.json results
            rpsbml = rpSBML(path=glob.glob(os.path.join(tmp_output_folder, '*'))[0])
            json.dump({'asDict': rpsbml.asDict()}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['asDict'], self.data['rp2ToSBML_asDict'])


    def test_readRp2PathsCompounds(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rp2pathscompounds = self.rpreader.readRp2PathsCompounds(os.path.join('data', 'rpreader', 'rp2paths_compounds.tsv'))
            json.dump({'readrp2pathscompounds': rp2pathscompounds}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['readrp2pathscompounds'], self.data['readrp2pathscompounds'])

    def test_readRpPathways(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rp_transformations, rp_compounds = self.rpreader.readRpPathways(os.path.join('data', 'rpreader', 'rp_pathways.csv'))
            json.dump({'rp_transformations': rp_transformations, 'rp_compounds': rp_compounds}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['rp_transformations'], self.data['rp_transformations'])
            self.assertCountEqual(new_data['rp_compounds'], self.data['rp_compounds'])

    def test_readRp2PathsPathways(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rp2pathspathways = self.rpreader.readRp2PathsPathways(os.path.join('data', 'rpreader', 'rp2paths_pathways.csv'))
            json.dump({'readrp2pathspathways': rp2pathspathways}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['readrp2pathspathways'], self.data['readrp2pathspathways']) 

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
