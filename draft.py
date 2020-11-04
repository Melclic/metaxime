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

