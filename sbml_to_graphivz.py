import libsbml
from graphviz import Digraph

document = libsbml.readSBMLFromFile('/Users/melchior/Downloads/rpglobalscore/rp_4_5.sbml.xml')
model = document.getModel()

f = Digraph('rp_1_1', filename='rp_1_1.gv')

#add the reactions
f.attr('node', shape='box', rank='same')
for reaction in model.getListOfReactions():
	print('Making reaction: '+str(reaction.getId()))
	f.node(reaction.getId(), label=reaction.getId())

f.attr('node', shape='ellipse')
for species in model.getListOfSpecies():
	print('Making species: '+str(species.getId()))
	f.node(species.getId(), label=species.getName())

for reaction in model.getListOfReactions():
	for reac in reaction.getListOfReactants():
		print('Making edge: '+str(reac.species)+' --> '+str(reaction.getId()))
		f.edge(reac.species, reaction.getId())
	for prod in reaction.getListOfProducts():
		f.edge(reaction.getId(), prod.species)
		print('Making edge: '+str(reaction.getId())+' --> '+str(prod.species))

f.view()






f.attr(rankdir='LR', size='8,5')
f.attr('node', shape='doublecircle')
f.node('LR_0')
f.node('LR_3')
f.node('LR_4')
f.node('LR_8')

f.attr('node', shape='circle')
f.edge('LR_0', 'LR_2', label='SS(B)')
f.edge('LR_0', 'LR_1', label='SS(S)')
f.edge('LR_1', 'LR_3', label='S($end)')
f.edge('LR_2', 'LR_6', label='SS(b)')
f.edge('LR_2', 'LR_5', label='SS(a)')
f.edge('LR_2', 'LR_4', label='S(A)')
f.edge('LR_5', 'LR_7', label='S(b)')
f.edge('LR_5', 'LR_5', label='S(a)')
f.edge('LR_6', 'LR_6', label='S(b)')
f.edge('LR_6', 'LR_5', label='S(a)')
f.edge('LR_7', 'LR_8', label='S(b)')
f.edge('LR_7', 'LR_5', label='S(a)')
f.edge('LR_8', 'LR_6', label='S(b)')
f.edge('LR_8', 'LR_5', label='S(a)')

f.view()



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
	for reac in reaction.getListOfReactants():
		g.add_edge(reac.species, reaction.getId())
	for prod in reaction.getListOfProducts():
		g.add_edge(reaction.getId(), prod.species)

nx.drawing.nx_pydot.write_dot(g, '/Users/melchior/Downloads/test_draw')