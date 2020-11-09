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

    @classmethod
    def setUpClass(cls):
        #load a rpSBML file
        cls.rpsele = rpSelenzyme(path=os.path.join('data', 'rpselenzyme', 'rpsbml.xml'))
        cls.data = json.load(open(os.path.join('data', 'rpselenzyme', 'data.json'), 'r'))

    def test_singleReactionRule(self):
        reac = self.rpsele.model.getReaction('RP1')
        annot = self.rpsele.readBRSYNTHAnnotation(a.getAnnotation())
        self.assertDictEqual(self.rpsele.singleReactionRule(annot['smiles'], 83333), self.data['singlereactionrule'])

    def test_run(self):
        self.assertDictEqual(self.rpsele.run(83333), self.data['run'])
        #add test to compare files hash when adding the results to rpSBML

if __name__ == '__main__':
    unittest.main()
