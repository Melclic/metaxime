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
        self.assertCountEqual(self.rpgraph.onlyConsumedSpecies(True, True), ['MNXM100__64__MNXC3'])
        self.assertCountEqual(self.rpgraph.onlyConsumedSpecies(True, False), ['MNXM100__64__MNXC3'])

    #onlyProducedSpecies
    def test_onlyProducedSpecies(self):
        self.assertCountEqual(self.rpgraph.onlyProducedSpecies(True, True), ['TARGET_0000000001__64__MNXC3'])
        self.assertCountEqual(self.rpgraph.onlyProducedSpecies(True, False), ['TARGET_0000000001__64__MNXC3'])

    #compare
    def test_compare(self):
        rpgraph = rpGraph(model_name='test', path=os.path.join('data', 'rpgraph', 'rpsbml_compare.xml'))
        self.assertAlmostEqual(rpGraph.compare(self.rpgraph, rpgraph), 0.23529411764705888)

    #orderedRetroReactions

    #asDict
    def test_networkDict(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            net_dict = self.rpgraph.networkDict()
            json.dump({'asdict': net_dict}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['asdict'], self.data['asdict'])

    #orderedRetroReactions


if __name__ == '__main__':
    unittest.main()
