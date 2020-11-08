import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpSelenzyme

class TestRPSelenzyme(unittest.TestCase):
    def setUp(self):
        #load a rpSBML file
        self.rpsele = rpSelenzyme(path=os.path.join('data', 'rpselenzyme', 'rpsbml.xml'))
        self.data = json.load(open(os.path.join('data', 'rpselenzyme', 'data.json'), 'r'))

    def test_singleReactionRule(self):
        reac = self.rpsele.model.getReaction('RP1')
        annot = self.rpsele.readBRSYNTHAnnotation(a.getAnnotation())
        self.assertDictEqual(self.rpsele.singleReactionRule(annot['smiles'], 83333), self.data['singlereactionrule'])

    def test_run(self):
        self.assertDictEqual(self.rpsele.run(83333), self.data['run'])

if __name__ == '__main__':
    unittest.main()
