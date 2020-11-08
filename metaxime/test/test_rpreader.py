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
    def setUp(self):
        #self.data = json.load(open(os.path.join('data', 'rpreader', 'rp2paths_compounds.csv'), 'r'))
        #rp2paths_compounds.csv rp2paths_pathways.csv  rp_pathways.csv

    def test_rp2ToSBML(self):
        rpcache = rpCache()
        rpcache.populateCache()
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            os.mkdir(os.path.join(tmp_output_folder, 'output'))
            rpreader = rpReader(rpcache)

    def readRp2Compounds

    def readRpPathways

    def readRp2PathsPathways

    def _gatherSpeciesInfo

    def _completeReac

    def _addCofactorSpecies

    def _updateStoichio

    def _addCofactorsStep

    def _addCofactors

    def _parseTSV

    def TSVtoSBML
