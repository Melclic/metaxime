from metaxime import rpCache
import os
rpcache = rpCache()
from metaxime import rpReader
rpReader.rp2ToCollection(os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpreader', 'rp_pathways.csv'),
                         os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpreader', 'rp2paths_compounds.tsv'),
                         os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpreader', 'rp2paths_pathways.csv'),
                         os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpreader', 'test.rpcol'),
                         rpcache=rpcache)

from metaxime import rpFBA
rpFBA.runCollection(os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpreader', 'test.rpcol'),
                    os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpfba', 'gem.xml'),
                    os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpfba', 'test.rpcol'),
                    num_workers=1,
                    rpcache=rpcache )

from metaxime import rpEquilibrator

rpEquilibrator.runCollection(os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpfba', 'test.rpcol'),
                    		 os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpequilibrator', 'test.rpcol'),
                    		 rpcache=rpcache)
from metaxime import rpSelenzyme
rpSelenzyme.runCollection(os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpequilibrator', 'test.rpcol'),
						  83333,
						  os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpselenzyme', 'test.rpcol'),
						  rpcache=rpcache)

from metaxime import rpGlobalScore
rpGlobalScore.runCollection(os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpselenzyme', 'test.rpcol'),
							os.path.join('/home/mdulac/workspace/melclic/metaxime/metaxime/test', 'data', 'rpglobalscore', 'test.rpcol'),
							rpcache=rpcache)

from metaxime import rpMerge
rpMerge.mergeSBMLFiles('/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpfba/rpsbml_collection/models/rp_129_6_rpsbml.xml',
					   '/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpfba/gem.xml',
					   '/home/mdulac/workspace/melclic/metaxime/metaxime/test/data/rpfba/test.sbml')
