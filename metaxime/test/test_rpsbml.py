import unittest
import os
import json
import hashlib
import tempfile
import libsbml
import sys

sys.path.insert(0, '../..')

from metaxime import rpSBML

class TestRPSBML(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        #load a rpSBML file
        self.rpsbml = rpSBML('test', path=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'rpsbml', 'rpsbml.xml'))
        self.gem = rpSBML('gem', path=os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'rpsbml', 'gem.xml'))
        self.data = json.load(open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', 'rpsbml', 'data.json'), 'r'))
        self.maxDiff = None #to be able to compare large dict

    def test_isRPsbml(self):
        self.assertTrue(self.rpsbml._isRPsbml())
        self.assertFalse(self.gem._isRPsbml())

    def test_getGroupsMembers(self):
        self.assertCountEqual(self.rpsbml.getGroupsMembers('rp_pathway'), ['RP1'])
        self.assertRaises(AttributeError, self.gem.getGroupsMembers, 'rp_pathway')

    def test_computeMeanRulesScore(self):
        self.assertAlmostEqual(self.rpsbml._computeMeanRulesScore(), 0.6702613998846042)
        self.assertRaises(AttributeError, self.gem._computeMeanRulesScore)

    def test_dictRPpathway(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            json.dump({'rppathway': self.rpsbml._dictRPpathway()}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['rppathway'], self.data['dictrppathway'])
            self.assertRaises(AttributeError, self.gem._dictRPpathway)

    def test_nameToSbmlId(self):
        self.assertEqual(self.rpsbml._nameToSbmlId('test123-_!"£$%^&*(){}@~><>?'), 'test123___________________')

    def test_genMetaID(self):
        self.assertEqual(self.rpsbml._genMetaID('test123'), 'cc03e747a6afbbcbf8be7668acfebee5')

    def test_readRPrules(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            json.dump({'rprules': self.rpsbml.readRPrules()}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['rprules'], self.data['readrprules'])
            self.assertRaises(AttributeError, self.gem.readRPrules)

    def test_readRPspecies(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            json.dump({'rpspecies': self.rpsbml.readRPspecies()}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['rpspecies'], self.data['readrpspecies'])
            self.assertRaises(AttributeError, self.gem.readRPspecies)

    def test_readUniqueRPspecies(self):
        self.assertCountEqual(self.rpsbml.readUniqueRPspecies(), ['TARGET_0000000001__64__MNXC3', 'MNXM100__64__MNXC3'])
        self.assertRaises(AttributeError, self.gem.readUniqueRPspecies)

    def test_readTaxonAnnotation(self):
        self.assertCountEqual(self.gem.readTaxonomy(), ['511145'])
        self.assertCountEqual(self.rpsbml.readTaxonomy(), [])

    def test_defaultBothAnnot(self):
        self.assertEqual(self.rpsbml._defaultBothAnnot('test'), self.data['defaultbothannot'])

    def test_defaultBRSynthAnnot(self):
        self.assertEqual(self.rpsbml._defaultBRSynthAnnot('test'), self.data['defaultbrsynthannot'])

    def test_defaultMIRIAMAnnot(self):
        self.assertEqual(self.rpsbml._defaultMIRIAMAnnot('test'), self.data['defaultmiriamannot'])

    def test_updateBRSynthPathway(self):
        rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
        tmp_dict = rpsbml.asDict()
        tmp_dict['pathway']['brsynth']['test'] = {'value': 99}
        rpsbml.updateBRSynthPathway(tmp_dict)
        self.assertEqual(rpsbml.asDict()['pathway']['brsynth']['test']['value'], 99)

    def test_addMIRIAMinchiKey(self):
        #no inchikeys should be added
        rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
        self.assertTrue(rpsbml.addMIRIAMinchiKey())
        self.assertEqual(rpsbml.readMIRIAMAnnotation(rpsbml.model.getSpecies('MNXM100__64__MNXC3').getAnnotation())['inchikey'][0], 'GVVPGTZRZFNKDS-JXMROGBWSA-K')
        #inchikeys are added
        gem = rpSBML('test', path=os.path.join('data', 'rpsbml', 'gem.xml'))
        self.assertTrue(gem.addMIRIAMinchiKey())
        self.assertEqual(gem.readMIRIAMAnnotation(gem.model.getSpecies('M_2pg_c').getAnnotation())['inchikey'][0], 'GXIURPTVHJPJLF-UWTATZPHSA-K') 

    #def test_overwriteRPannot(self):

    #def addUpdateBRSynth(self):

    #def addUpdateMIRIAM(self):

    def test_asDict(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            json.dump({'asdict': self.rpsbml.asDict()}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['asdict'], self.data['asdict'])
            self.assertRaises(AttributeError, self.gem.asDict)

    #def test_readSBML(self):

    #def test_writeSBML(self):

    def test_findCreateObjective(self):
        #the find part
        self.assertEqual(self.rpsbml.findCreateObjective(['RP1_sink'], [1.0]), 'obj_RP1_sink')
        #the create part
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
            self.assertEqual(rpsbml.findCreateObjective(['RP2'], [1.0]), 'obj_RP2')


    def test_readMIRIAMAnnotation(self):
        self.assertDictEqual(self.rpsbml.readMIRIAMAnnotation(self.rpsbml.model.getReaction('RP1').getAnnotation()), {'ec-code': ['4.2.3.16', '4.2.3.20']})
        self.assertDictEqual(self.gem.readMIRIAMAnnotation(self.gem.model.getReaction('R_ALATA_D2').getAnnotation()), {'bigg': ['ALATA_D2'], 'biocyc': ['RXN0-5240'], 'kegg': ['R01147'], 'metanetx': ['MNXR95697'], 'rhea': ['28562', '28563', '28564', '28565']})

    def test_readBRSYNTHAnnotation(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            json.dump({'readbrsynthannotation': self.rpsbml.readBRSYNTHAnnotation(self.rpsbml.model.getReaction('RP1').getAnnotation())},
                    open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['readbrsynthannotation'], self.data['readbrsynthannotation'])

    def test_readReactionSpecies(self):
        self.assertDictEqual(self.rpsbml.readReactionSpecies(self.rpsbml.model.getReaction('RP1')), {'left': {'MNXM100__64__MNXC3': 1}, 'right': {'TARGET_0000000001__64__MNXC3': 1}})

    def test_speciesExists(self):
        self.assertTrue(self.rpsbml.speciesExists('MNXM100'))
        self.assertFalse(self.rpsbml.speciesExists('test'))

    def test_isSpeciesProduct(self):
        self.assertTrue(self.rpsbml.isSpeciesProduct('TARGET_0000000001__64__MNXC3'))
        self.assertFalse(self.rpsbml.isSpeciesProduct('MNXM100__64__MNXC3'))

    def test_outPathsDict(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            json.dump({'outpathsdict': self.rpsbml.outPathsDict()}, open(os.path.join(tmp_output_folder, 'data.json'), 'w'))
            new_data = json.load(open(os.path.join(tmp_output_folder, 'data.json')))
            self.assertDictEqual(new_data['outpathsdict'], self.data['outpathsdict'])
            self.assertRaises(AttributeError, self.gem.outPathsDict)

    def test_compareBRSYNTHAnnotations(self):
        self.assertTrue(self.rpsbml.compareBRSYNTHAnnotations(self.rpsbml.model.getSpecies('MNXM100__64__MNXC3').getAnnotation(), self.rpsbml.model.getSpecies('MNXM100__64__MNXC3').getAnnotation()))
        self.assertFalse(self.rpsbml.compareBRSYNTHAnnotations(self.rpsbml.model.getSpecies('MNXM100__64__MNXC3').getAnnotation(), self.rpsbml.model.getSpecies('TARGET_0000000001__64__MNXC3').getAnnotation()))
        self.assertFalse(self.gem.compareBRSYNTHAnnotations(self.gem.model.getSpecies('M_2pg_c').getAnnotation(), self.gem.model.getSpecies('M_13dpg_c').getAnnotation()))

    def test_compareMIRIAMAnnotations(self):
        self.assertTrue(self.rpsbml.compareMIRIAMAnnotations(self.rpsbml.model.getSpecies('MNXM100__64__MNXC3').getAnnotation(), self.rpsbml.model.getSpecies('MNXM100__64__MNXC3').getAnnotation()))
        self.assertFalse(self.rpsbml.compareMIRIAMAnnotations(self.rpsbml.model.getSpecies('MNXM100__64__MNXC3').getAnnotation(), self.rpsbml.model.getSpecies('TARGET_0000000001__64__MNXC3').getAnnotation()))
        
    def test_genericModel(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            new = rpSBML('test')
            new.genericModel('test_name', 'test_id', 'MNXC3')
            new.writeSBML(os.path.join(tmp_output_folder, 'test.xml'))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), 'f1de251a278b40ef0448579e308edf9d')

    def test_createModel(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            new = rpSBML('test')
            new.createModel('test_name', 'test_id')
            new.writeSBML(os.path.join(tmp_output_folder, 'test.xml'))           
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), '79b9063da98b5bd8286afab444f52227')
            
    def test_createCompartment(self):
        #test with no comp xref
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            new = rpSBML('test')
            new.createModel('test_name', 'test_id')
            new.createCompartment('test')
            new.writeSBML(os.path.join(tmp_output_folder, 'test.xml'))           
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), '87fd0a156f470799bdcf6295decb5faa')
        #test with xref
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            new = rpSBML('test')
            new.createModel('test_name', 'test_id')
            new.createCompartment('MNXC3')
            new.writeSBML(os.path.join(tmp_output_folder, 'test.xml'))           
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), 'c27be1b46e4c0f10aaf86362364749e7')
        
    def test_createUnitDefinition(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            new = rpSBML('test')
            new.createModel('test_name', 'test_id')
            unit_def = new.createUnitDefinition('mmol_per_gDW_per_hr')
            new.createUnit(unit_def, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
            new.createUnit(unit_def, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
            new.createUnit(unit_def, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
            new.writeSBML(os.path.join(tmp_output_folder, 'test.xml'))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), '7e98e6667c2de2fa9078795121203c99')
        
    def test_createReturnFluxParameter(self):
        #create feature
        new = rpSBML('test')
        new.createModel('test_name', 'test_id')
        param = new.createReturnFluxParameter(8888.0)
        self.assertEqual(param.id, 'B_8888_0')
        self.assertEqual(param.value, 8888.0)

    def test_createReaction(self):
        #TODO: add detection
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
            step = {'rule_id': None,
                   'left': {'MNXM89557': 1},
                   'right': {},
                   'step': None,
                   'sub_step': None,
                   'path_id': None,
                   'transformation_id': None,
                   'rule_score': None,
                   'rule_ori_reac': None}
            rpsbml.createReaction('test',
                                  999999.0,
                                  0.0,
                                  step,
                                  'MNXC3',
                                  xref={'ec': ['1.1.1.1']})  
            rpsbml.writeSBML(os.path.join(tmp_output_folder, 'test.xml')) 
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), '6c832eb02bb3f1849d2e0b754e44e748')

    def test_createSpecies(self):
        #create new species
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
            spe = rpsbml.createSpecies('test', 'MNXC3')
            self.assertEqual(spe.getId(), 'test')
            rpsbml.writeSBML(os.path.join(tmp_output_folder, 'test.xml'))            
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), '8222edc93f5c1430231d7e16afbf3079')
        #recover already existing species
        rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
        spe = rpsbml.createSpecies('MNXM6', 'MNXC3')
        self.assertEqual(spe.getId(), 'MNXM6')
            
    def test_createGroup(self):
        #create new group
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
            gro = rpsbml.createGroup('test')
            self.assertEqual(gro.getId(), 'test')
            rpsbml.writeSBML(os.path.join(tmp_output_folder, 'test.xml'))            
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), '19c61349df77d97118210a359710dbbf')
        #recover already existing group
        rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
        gro = rpsbml.createGroup('rp_pathway', 'MNXC3')
        self.assertEqual(gro.getId(), 'rp_pathway')
        
    #def test_createGene(self):

    def test_createMultiFluxObj(self):
        #create new group
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
            flux_obj = rpsbml.createMultiFluxObj('test', ['RP1'], [1.0])
            self.assertEqual(flux_obj.getId(), 'test')
            rpsbml.writeSBML(os.path.join(tmp_output_folder, 'test.xml'))
            self.assertEqual(hashlib.md5(open(os.path.join(tmp_output_folder, 'test.xml'), 'rb').read()).hexdigest(), 'c4647ea976dbf1eebbf0e7400af4b65d')
        #recover already existing group
        rpsbml = rpSBML('test', path=os.path.join('data', 'rpsbml', 'rpsbml.xml'))
        flux_obj = rpsbml.createMultiFluxObj('obj_fraction', ['RP1'], [1.0])
        self.assertEqual(flux_obj.getId(), 'obj_fraction')
        

if __name__ == '__main__':
    unittest.main()
