import unittest
import sys

sys.path.insert(0, '../..')

from metaxime import rpReader

class TestRPSBML(unittest.TestCase):
    def setUp(self):
        self.data = json.load(open(os.path.join('data', 'rpsbml', 'data.json'), 'r'))

    def test_isRPsbml(self):
