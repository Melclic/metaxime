import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpExtractSink

class TestRPSelenzyme(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        #load a rpSBML file
        #self.data = json.load(open(os.path.join('data', 'rpselenzyme', 'data.json'), 'r'))
        self.maxDiff = None

    #def _convertToCobra

    #def _reduce_model

    #def _removeDeadEnd

    def test_genSink(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpExtractSink.genSink(os.path.join('data', 'rpselenzyme', 'data.json'), os.path.join(tmp_output_folder, 'tmp.csv'), remove_dead_end=True)

if __name__ == '__main__':
    unittest.main()
