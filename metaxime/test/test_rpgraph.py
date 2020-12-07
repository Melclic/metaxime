import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpGraph

class TestRPGraph(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        #load a rpSBML file
        self.rpgraph = rpGraph(model_name='test', path=os.path.join('data', 'rpgraph', 'rpsbml.xml'))
        self.data = json.load(open(os.path.join('data', 'rpgraph', 'data.json'), 'r'))
        self.maxDiff = None


    #_makeCompareGraphs

    #_makeGraph

    #_recursiveReacSuccessors

    #_recursiveReacPredecessors

    #_recursiveReacPredecessors

    #onlyConsumedSpecies
    def test_onlyConsumedSpecies(self):
        self.assertCountEqual(self.rpgraph.onlyConsumedSpecies(True, True), ['MNXM89557__64__MNXC3', 'MNXM1__64__MNXC3', 'MNXM6__64__MNXC3', 'MNXM3__64__MNXC3'])
        self.assertCountEqual(self.rpgraph.onlyConsumedSpecies(True, False), ['MNXM89557__64__MNXC3', 'MNXM1__64__MNXC3'])

    #onlyProducedSpecies
    def test_onlyProducedSpecies(self):
        self.assertCountEqual(self.rpgraph.onlyProducedSpecies(True, True), ['TARGET_0000000001__64__MNXC3', 'MNXM9__64__MNXC3', 'MNXM5__64__MNXC3', 'MNXM7__64__MNXC3', 'MNXM20__64__MNXC3', 'MNXM13__64__MNXC3'])
        self.assertCountEqual(self.rpgraph.onlyProducedSpecies(True, False), ['TARGET_0000000001__64__MNXC3'])

    #compare
    def test_compare(self):
        rpgraph = rpGraph(model_name='test', path=os.path.join('data', 'rpgraph', 'rpsbml_compare.xml'))
        self.assertEqual(rpGraph.compare(self.rpgraph, rpgraph), 0.75)

    #orderedRetroReactions

    #exportJSON
    def test_exportJSON(self):
        self.assertDictEqual(self.rpgraph.exportJSON(), self.data['exportjson'])

    #orderedRetroReactions


if __name__ == '__main__':
    unittest.main()
