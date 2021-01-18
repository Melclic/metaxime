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
    def setUpClass(self):
        #load a rpSBML file
        self.rpsele = rpSelenzyme(path=os.path.join('data', 'rpselenzyme', 'rpsbml.xml'))
        #self.data = json.load(open(os.path.join('data', 'rpselenzyme', 'data.json'), 'r'))
        self.maxDiff = None

    def test_singleReactionRule(self):
        reac = self.rpsele.model.getReaction('RP1')
        annot = self.rpsele.readBRSYNTHAnnotation(reac.getAnnotation())
        uniprot = self.rpsele.singleReactionRule(annot['smiles'], 83333)
        self.assertEqual(uniprot['P07773'], 99.9)
        #self.assertDictEqual(self.rpsele.singleReactionRule(annot['smiles'], 83333), self.data['singlereactionrule'])

    '''
    def test_run(self):
        self.assertDictEqual(self.rpsele.run(83333), self.data['run'])
        #add test to compare files hash when adding the results to rpSBML
    '''

if __name__ == '__main__':
    unittest.main()
