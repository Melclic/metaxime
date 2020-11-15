import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpGlobalScore

class TestRPSBML(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        #load a rpSBML file
        cls.rpglobalscore = rpGlobalScore(path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))

    def test_calculateGlobalScore(self):
        score = self.rpglobalscore.calculateGlobalScore()
        self.assertEqual(score, 0.0)

if __name__ == '__main__':
    unittest.main()
