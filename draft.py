########################## MAC ###########################

from metaxime import rpSBML
import hashlib
import json
import libsbml

rpsbml = rpSBML(model_name='test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpsbml/rpsbml.xml')
gem = rpSBML(model_name='test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpsbml/gem.xml')

rpsbml1 = rpSBML(model_name='test', path='/Users/melchior/Downloads/rpglobalscore_77/rp_9_1.sbml.xml')
#rpsbml1.updateBRSynthPathway(rpsbml.asDict())
#rpsbml1.writeSBML('/Users/melchior/Downloads/test.xml')

rpsbml1_dict = rpsbml1.asDict()
rpsbml1_dict['pathway']['brsynth']['global_score']['value']

rpsbml_dict = rpsbml.asDict()
rpsbml_dict['pathway']['brsynth']['global_score']['value']

rpsbml1.updateBRSynthPathway(rpsbml.asDict())

rpsbml1_mod_dict = rpsbml1.asDict()
rpsbml1_dict['pathway']['brsynth']['global_score']['value']
rpsbml_dict['pathway']['brsynth']['global_score']['value']
rpsbml1_mod_dict['pathway']['brsynth']['global_score']['value']

'/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpgraph/data.json'

from metaxime import rpGraph
import hashlib
import json
import libsbml

rpgraph = rpGraph(model_name='test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpgraph/rpsbml.xml')
rpgraph_compare = rpGraph(model_name='test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpgraph/rpsbml_compare.xml')
rpGraph.compare(rpgraph, rpgraph_compare)


rpsbml = rpSBML(model_name='test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpmerge/rpsbml.xml')
gem = rpSBML(model_name='test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpmerge/gem.xml')


from metaxime import rpMerge
import hashlib
rpMerge.mergeSBMLFiles('/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpmerge/rpsbml.xml', '/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpmerge/gem.xml', '/Users/melchior/Downloads/test.xml')

from metaxime import rpFBA
import hashlib
rpfba = rpFBA(rpsbml_path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpfba/rpsbml.xml', gem_sbml_path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpfba/gem.xml')

rpfba.writeSBML('/Users/melchior/Downloads/merged.xml')
hashlib.md5(open('/Users/melchior/Downloads/merged.xml', 'rb').read()).hexdigest()
'0e755a7ae4605279df728b5dab176181'


from metaxime import rpFBA
import hashlib
rpfba = rpFBA(rpsbml_path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpfba/rpsbml.xml', gem_sbml_path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpfba/gem.xml')
rpfba.runFractionReaction('biomass', 1.0, 'RP1_sink', 1.0)
rpfba.writeSBML('/Users/melchior/Downloads/merged.xml', False)


from metaxime import rpEquilibrator


rpgraph.compare(rpgraph_compare)


data = json.load(open('/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpgraph/data.json', 'r'))


json.dump(data, open('/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpgraph/data.json', 'w'))







data = json.load(open('/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpsbml/data.json', 'r'))








json.dump(data, open('/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/data.json', 'w'))




from metaxime import rpDraw
rpdraw = rpDraw(path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpsbml.xml')
svg, resG, mod_pos, reac_cofactors_name, nodes_attach_locs = rpdraw.drawSVG('/Users/melchior/Downloads/test.svg')



import rpSBML
import rpGraph
import rpDraw

rpsbml = rpSBML.rpSBML('test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpsbml.xml')
rpgraph = rpGraph.rpGraph(rpsbml)
rpdraw = rpDraw.rpDraw()

svg, resG, mod_pos, reac_cofactors_name, nodes_attach_locs = rpdraw.drawsvg(rpgraph, 'TARGET_0000000001__64__MNXC3')

open('/Users/melchior/Downloads/test.svg', 'w').write(svg)


###################### LINUX ##################################

from metaxime import rpReader
rpReader.rp2ToCollection('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpreader/rp_pathways.csv', '/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpreader/rp2paths_compounds.csv', '/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpreader/rp2paths_pathways.csv', '/home/mdulac/Downloads/test_out.rpcol')


from metaxime import rpEquilibrator
rpeq = rpEquilibrator(path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpequilibrator/rpsbml.xml')
spe = rpeq.model.getSpecies('MNXM89557__64__MNXC3')
rpeq._makeSpeciesStr(spe)




json.dump(data, open('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpreader/data.json', 'w'))

from metaxime import rpSBML
import hashlib
import json
import libsbml

rpsbml = rpSBML('test', path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpsbml.xml')
gem = rpSBML('gem', path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/gem.xml')

data = json.load(open('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/data.json', 'r'))




data = json.load(open('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpselenzyme/data.json', 'r'))


from metaxime import rpSelenzyme
#rpsele = rpSelenzyme('/home/mdulac/workspace/melclic/metaxime/metaxime/input_cache/rpselenzyme_data.tar.xz', path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpsbml/rpsbml.xml')
rpsele = rpSelenzyme(path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpsbml/rpsbml.xml')
res = rpsele.run(83333)


from metaxime import rpSelenzyme
rpsele = rpSelenzyme(path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpsbml/rpsbml.xml')
a = rpsele.model.getReaction('RP1')
b = rpsele.readBRSYNTHAnnotation(a.getAnnotation())
uniprotID_score = rpsele.singleReactionRule(b['smiles'], 83333)

json.dump(data, open('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpselenzyme/data.json', 'w'))


uniprot_aaLenght = {}
with open('sel_len.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    next(csv_reader)
    for row in csv_reader:
        uniprot_aaLenght[row[0].split('|')[1]] = int(row[1])


new = rpSBML('test')
new.genericModel('test_name', 'test_id', 'MNXC3')
new.writeSBML('/home/mdulac/Downloads/genericmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/genericmodel.xml', 'rb').read()).hexdigest()


new = rpSBML('test')
new.createModel('test_name', 'test_id')
new.writeSBML('/home/mdulac/Downloads/createmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/createmodel.xml', 'rb').read()).hexdigest()


new = rpSBML('test')
new.createModel('test_name', 'test_id')
new.createCompartment('test')
new.writeSBML('/home/mdulac/Downloads/createmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/createmodel.xml', 'rb').read()).hexdigest()
87fd0a156f470799bdcf6295decb5faa

new = rpSBML('test')
new.createModel('test_name', 'test_id')
new.createCompartment('MNXC3')
new.writeSBML('/home/mdulac/Downloads/createmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/createmodel.xml', 'rb').read()).hexdigest()
c27be1b46e4c0f10aaf86362364749e7


new = rpSBML('test')
new.createModel('test_name', 'test_id')
unit_def = new.createUnitDefinition('mmol_per_gDW_per_hr')
new.createUnit(unit_def, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
new.createUnit(unit_def, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
new.createUnit(unit_def, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
new.writeSBML('/home/mdulac/Downloads/createmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/createmodel.xml', 'rb').read()).hexdigest()


new = rpSBML('test')
new.createModel('test_name', 'test_id')
unit_def = new.createUnitDefinition('mmol_per_gDW_per_hr')
new.createUnit(unit_def, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
new.createUnit(unit_def, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
new.createUnit(unit_def, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
new.writeSBML('/home/mdulac/Downloads/createmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/createmodel.xml', 'rb').read()).hexdigest()

param = rpsbml.createReturnFluxParameter(None, parameter_id='B_999999')
rpsbml.createReturnFluxParameter(999999)
new = rpSBML('test')
new.createModel('test_name', 'test_id')
param = new.createReturnFluxParameter(8888)
param.id = 'B_8888'
param.value = 8888.0

new = rpSBML('test')
new.createModel('test_name', 'test_id')
unit_def = new.createUnitDefinition('mmol_per_gDW_per_hr')
new.createUnit(unit_def, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
new.createUnit(unit_def, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
new.createUnit(unit_def, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
new.writeSBML('/home/mdulac/Downloads/createmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/createmodel.xml', 'rb').read()).hexdigest()

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
rpsbml.writeSBML('/home/mdulac/Downloads/createmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/createmodel.xml', 'rb').read()).hexdigest()
4d73b3390ec35368d4fbc63b961908d1


rpsbml.createSpecies('test', 'MNXC3')
rpsbml.writeSBML('/home/mdulac/Downloads/createmodel.xml')
hashlib.md5(open('/home/mdulac/Downloads/createmodel.xml', 'rb').read()).hexdigest()

rpsbml.createSpecies('MNXM6', 'MNXC3')


rpsbml.createMultiFluxObj('test', ['RP1'], [1.0])


fbc_plugin = rpsbml.model.getPlugin('fbc')





from metaxime import rpSBML
import hashlib
import json
import libsbml

#data = json.load(open('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/data.json', 'r'))
data = {}
rpsbml = rpSBML('test', path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpsbml.xml')

data['dictrppathway'] = rpsbml._dictRPpathway()
data['readrprules'] = rpsbml.readRPrules()
data['readrpspecies'] = rpsbml.readRPspecies()
data['defaultbothannot'] = rpsbml._defaultBothAnnot('test')
data['defaultbrsynthannot'] = rpsbml._defaultBRSynthAnnot('test')
data['defaultmiriamannot'] = rpsbml._defaultMIRIAMAnnot('test')
data['asdict'] = rpsbml.asDict()
data['readbrsynthannotation'] = rpsbml.readBRSYNTHAnnotation(rpsbml.model.getReaction('RP1').getAnnotation())


json.dump(data, open('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/data.json', 'w'))











from metaxime import rpDraw
rpdraw = rpDraw(path='')
svg, resG, mod_pos, reac_cofactors_name, nodes_attach_locs = rpdraw.drawSVG('/home/mdulac/Downloads/test.svg')



import rpSBML
import rpGraph
import rpDraw

rpsbml = rpSBML.rpSBML('test', path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpsbml.xml')
rpgraph = rpGraph.rpGraph(rpsbml)
rpdraw = rpDraw.rpDraw()

svg, resG, mod_pos, reac_cofactors_name, nodes_attach_locs = rpdraw.drawsvg(rpgraph, 'TARGET_0000000001__64__MNXC3')

open('/home/mdulac/Downloads/test.svg', 'w').write(svg)

