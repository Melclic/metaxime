########################## MAC ###########################

from metaxime import rpSBML
import json

rpsbml = rpSBML('test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpsbml.xml')
gem = rpSBML('test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/gem.xml')

data = json.load(open('/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/data.json', 'r'))

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


from metaxime import rpSBML
import hashlib
import json
import libsbml

rpsbml = rpSBML('test', path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpsbml.xml')
gem = rpSBML('gem', path='/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/gem.xml')

data = json.load(open('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/data.json', 'r'))


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

