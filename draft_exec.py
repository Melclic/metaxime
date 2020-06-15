import rpSBML
import numpy as np
import networkx as nx
import rpGraph

pathway_id = 'rp_pathway'

rpsbml = rpSBML.rpSBML('rp_4_5')
rpsbml.readSBML('/Users/melchior/Downloads/rpglobalscore/rp_4_5.sbml.xml')

rpgraph = rpGraph.rpGraph(rpsbml, pathway_id)

rpgraph.endSpecies()

rpgraph.startSpecies()

rpgraph.orderedRetroReaction()

rpgraph.orderedReaction()
