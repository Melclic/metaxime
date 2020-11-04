import unittest
import os
import json
import hashlib
import tempfile

from metaxime import rpSBML

class TestRPSBML(unittest.TestCase):
    def setUp(self):
        #load a rpSBML file
        self.rpsbml = rpSBML('test', path=os.path.join(os.path.dirname(__file__), 'data', 'rpsbml.xml'))
        self.gem = rpSBML('gem', path=os.path.join(os.path.dirname(__file__), 'data', 'gem.xml'))
        self.data = json.load(os.path.join(os.path.dirname(__file__), 'data', 'data.json'))

    def test_getGroupsMembers(self):
        self.assertCountEqual(self.rpsbml.getGroupsMembers('rp_pathway'), ['RP1', 'RP2', 'RP3'])
        self.asserRaises(AttributeError, self.gem.getGroupsMembers('rp_pathway'))

    def test_computeMeanRulesScore(self):
        self.assertAlmostEqual(self.rpsbml._computeMeanRulesScore(), 0.5684564101634014, places=7, message='Equal to 0.5684564101634014')
        self.asserRaises(AttributeError, self.gem._computeMeanRulesScore())

    def test_dictRPpathway(self):
        self.assertDictEqual(self.rpsbml._dictRPpathway(), self.data['dictrppathway'])
        self.asserRaises(AttributeError, self.gem._dictRPpathway())

    def test_nameToSbmlId(self):
        self.assertEqual(self.rpsbml._nameToSbmlId('test123-_!"£$%^&*(){}@~><>?'), 'test123___________________')

    def test_genMetaID(self):
        self.assertEqual(self.rpsbml._genMetaID('test123'), 'cc03e747a6afbbcbf8be7668acfebee5')

    def test_readRPrules(self):
        self.assertDictEqual(self.rpsbml.readRPrules(), self.data['readrprules'])
        self.asserRaises(AttributeError, self.gem.readRPrules())

    def test_readRPspecies(self):
        self.assertDictEqual(self.rpsbml.readRPspecies(), self.data['readrpspecies'])
        self.assertRaises(AttributeError, self.gem.readRPspecies())

    def test_readUniqueRPspecies(self):
        self.assertCountEqual(self.rpsbml.readUniqueRPspecies(), ['TARGET_0000000001__64__MNXC3', 'MNXM13__64__MNXC3', 'CMPD_0000000004__64__MNXC3', 'MNXM1__64__MNXC3', 'MNXM20__64__MNXC3', 'CMPD_0000000013__64__MNXC3', 'MNXM89557__64__MNXC3', 'MNXM5__64__MNXC3', 'MNXM7__64__MNXC3', 'MNXM9__64__MNXC3', 'MNXM6__64__MNXC3', 'MNXM3__64__MNXC3'])
        self.assertRaises(AttributeError, self.gem.readUniqueRPspecies())

    def test_readTaxonAnnotation(self):
        self.assertCountEqual(self.gem.readTaxonAnnotation(), ['511145'])
        self.assertIsNone(self.rpsbml.readTaxonAnnotation())

    def test_defaultBothAnnot(self):
        self.assertEqual(rpsbml._defaultBothAnnot('test'), self.data['defaultbothannot'])

    def test_defaultBRSynthAnnot(self):
        self.assertEqual(rpsbml._defaultBRSynthAnnot('test'), self.data['defaultbrsynthannot'])

    def test_defaultMIRIAMAnnot(self):
        self.assertEqual(rpsbml._defaultMIRIAMAnnot('test'), self.data['defaultmiriamannot'])

    def test_addMIRIAMinchiKey(self):
        #no inchikeys should be added
        self.assertTrue(self.rpsbml.addMIRIAMinchiKey())
        self.rpsbml.addMIRIAMinchiKey())
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            self.assertTrue(self.rpsbml.writeSBML(os.path.join(tmp_output_folder, 'test.sbml')))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.sbml'), 'rb').read()).hexdigest(), '6e9eafb99f4d2432188d6ecbda12a726')
        #inchikeys are added
        self.assertTrue(gem.addMIRIAMinchiKey())
        self.gem.addMIRIAMinchiKey())
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            self.assertTrue(self.gem.writeSBML(os.path.join(tmp_output_folder, 'test.sbml')))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.sbml'), 'rb').read()).hexdigest(), 'e7baa7d1a0a096a9dfe5c8160853df05')

    #def test_overwriteRPannot(self):

    #def addUpdateBRSynth(self):

    #def addUpdateMIRIAM(self):

    def test_asDict(self):
        self.assertDictEqual(self.rpsbml.asDict(), self.data['asdict'])
        self.asserRaises(AttributeError, self.gem._asDict())

    #def test_readSBML(self):

    #def test_writeSBML(self):

    def test_findCreateObjective(self):
        #the find part
        self.assertEqual(self.rpsbml.findCreateObjective(['RP1_sink'], [1.0]), 'obj_fraction')
        #the create part
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpsbml = rpSBML('test', path=os.path.join(os.path.dirname(__file__), 'data', 'rpsbml.xml'))
            self.assertEqual(rpsbml.findCreateObjective(['RP2'], [1.0]), 'obj_RP2')
            self.assertTrue(rpsbml.writeSBML(os.path.join(tmp_output_folder, 'test.sbml')))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.sbml'), 'rb').read()).hexdigest(), '66cb235c127e2bb07c6c13bea7bc6df2')

    def test_readMIRIAMAnnotation(self):
        self.assertDictEqual(self.rpsbml.readMIRIAMAnnotation(self.rpsbml.model.getReaction('RP1').getAnnotation()), {'ec-code': ['4.1.1.17']})
        self.assertDictEqual(self.gem.readMIRIAMAnnotation(self.gem.model.getReaction('R_ALATA_D2').getAnnotation()), {'bigg': ['ALATA_D2'], 'biocyc': ['RXN0-5240'], 'kegg': ['R01147'], 'metanetx': ['MNXR95697'], 'rhea': ['28562', '28563', '28564', '28565']})

    def test_readBRSYNTHAnnotation(self):
        self.assertDictEqual(self.rpsbml.readBRSYNTHAnnotation(self.rpsbml.model.getReaction('RP1').getAnnotation()), data['readbrsynthannotation'])
        self.assertDictEqual(self.gem.readBRSYNTHAnnotation(self.gem.model.getReaction('RP1').getAnnotation()), {})

    def test_readReactionSpecies(self):
        self.assertDictEqual(self.rpsbml.readReactionSpecies(self.rpsbml.model.getReaction('RP1')), {'left': {'CMPD_0000000004__64__MNXC3': 1, 'MNXM1__64__MNXC3': 1}, 'right': {'TARGET_0000000001__64__MNXC3': 1, 'MNXM13__64__MNXC3': 1}})

    def test_speciesExists(self):
        self.assertTrue(self.rpsbml.speciesExists('MNXM89557'))
        self.assertFalse(self.rpsbml.speciesExists('test'))

    def test_isSpeciesProduct(self):
        self.assertTrue(self.rpsbml.isSpeciesProduct('TARGET_0000000001__64__MNXC3'))
        self.assertFalse(self.rpsbml.isSpeciesProduct('MNXM1__64__MNXC3'))

    def test_outPathsDict(self):
        self.assertDictEqual(self.rpsbml.outPathsDict(), data['outpathsdict'])
        self.asserRaises(AttributeError, self.gem.outPathsDict())

    #def test_compareBRSYNTHAnnotations(self):

if __name__ == '__main__':
    unittest.main()
