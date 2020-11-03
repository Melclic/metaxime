import unittest
import os

from metaxime import rpSBML

class TestRPSBML(unittest.TestCase):
    def setUp(self):
        #load a rpSBML file
        self.rpsbml = rpSBML('test', path=os.path.join(os.path.dirname(__file__), 'data', 'test_rpsbml.xml'))
        self.dictrppathway = {'RP1': {'Reactants': {'AHLPHDHHMVZTML-UHFFFAOYSA-N': 1.0, 'GPRLSGONYQIRFK-UHFFFAOYSA-N': 1.0}, 'Products': {'KIDHWZJUCRJVML-UHFFFAOYSA-N': 1.0, 'CURLTUGMZLYLDI-UHFFFAOYSA-N': 1.0}}, 'RP2': {'Reactants': {'KABXUUFDPUOJMW-UHFFFAOYSA-N': 1.0, 'WHUUTDBJXJRKMK-VKHMYHEASA-M': 1.0}, 'Products': {'AHLPHDHHMVZTML-UHFFFAOYSA-N': 1.0, 'KPGXRSRHYNQIFN-UHFFFAOYSA-L': 1.0}}, 'RP3': {'Reactants': {'GPRLSGONYQIRFK-UHFFFAOYSA-N': 1.0, 'WHUUTDBJXJRKMK-VKHMYHEASA-M': 1.0, 'ACFIXJIJDZMPPO-NNYOXOHSSA-J': 1.0, 'ZKHQWZAMYRWXGA-KQYNXXCUSA-J': 1.0}, 'Products': {'KABXUUFDPUOJMW-UHFFFAOYSA-N': 1.0, 'XJLXINKUBYWONI-NNYOXOHSSA-K': 1.0, 'XTWYTFMLZFPYCI-KQYNXXCUSA-K': 1.0, 'NBIIXXVUZAFLBC-UHFFFAOYSA-L': 1.0}}}

    def test_computeMeanRulesScore(self):
        self.assertAlmostEqual(self.rpsbml._computeMeanRulesScore(), 0.5684564101634014, places=7, message='Equal to 0.5684564101634014')

    def test_dictRPpathway(self):
        self.assertDictEqual(rpsbml._dictRPpathway(), self.dictrppathway)

    def test_nameToSbmlId(self):
        self.assertEqual(rpsbml._nameToSbmlId('test123-_!"£$%^&*(){}@~><>?'), 'test123___________________')

    def test_genMetaID(self):
        self.assertEqual(rpsbml._genMetaID('test123'), 'cc03e747a6afbbcbf8be7668acfebee5')

if __name__ == '__main__':
    unittest.main()
