from metaxime import rpSBML
import json

rpsbml = rpSBML('test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/rpsbml.xml')
gem = rpSBML('test', path='/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/gem.xml')

data = json.load(open('/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/data.json', 'r'))

json.dump(data, open('/Users/melchior/workspace/melclic/metaxime/metaxime/test/data/data.json', 'w'))

