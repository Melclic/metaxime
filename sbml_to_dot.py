import libsbml
import networkx as nx

document = libsbml.readSBMLFromFile('/Users/melchior/Downloads/Galaxy85/rp_1_1.sbml.xml')
model = document.getModel()

g = nx.DiGraph()

def color_hex( color, prefix = '#' ):
	if len( color ) == 3:
		color = color + (255,)
	hexColor = prefix + ''.join( [ '%02x' % x for x in color ] )
	return hexColor

#add the species
color = (0,0,255,100)
for species in model.getListOfSpecies():
	g.add_node(species.getId(),
			   type = 'species', 
			   label = species.getName(),
			   color = 'red',
			   shape = 'circle',
			   fillcolor = color_hex(color),
			   style = 'filled',
			   fontsize = 8)

#add the reactions
for reaction in model.getListOfReactions():
	g.add_node(reaction.getId(),
			   type = 'reaction',
			   label = reaction.getId(),
			   shape = 'square')

for reaction in model.getListOfReactions():
	for reac in reaction.getListOfReactants():
		g.add_edge(reac.species, reaction.getId())
	for prod in reaction.getListOfProducts():
		g.add_edge(reaction.getId(), prod.species)

nx.drawing.nx_pydot.write_dot(g, '/Users/melchior/Downloads/test_draw')