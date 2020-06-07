import rpSBML
import numpy as np
import networkx as nx
import rpGraph

pathway_id = 'rp_pathway'

rpsbml = rpSBML.rpSBML('rp_4_5')
rpsbml.readSBML('/Users/melchior/Downloads/rpglobalscore/rp_4_5.sbml.xml')

G = nx.DiGraph()

species = [rpsbml.model.getSpecies(i) for i in rpsbml.readUniqueRPspecies(pathway_id)]
rp_dict = rpsbml.readRPspecies(pathway_id)
reactions = [rpsbml.model.getReaction(i) for i in list(rp_dict.keys())]


for spe in species:
	G.add_node(spe.getId(),
			   type='species')


for reac in reactions:
	G.add_node(reac.getId(),
			   type='reaction')


for reaction in reactions:
	for reac in reaction.getListOfReactants():
		G.add_edge(reac.species, reaction.getId())
	for prod in reaction.getListOfProducts():
		G.add_edge(reaction.getId(), prod.species)



#### species

for node_name in G.nodes.keys():
	node = G.nodes.get(node_name)
	if node['type']=='species':
		print('------------ '+str(node_name)+' ------------')
		print('successors:')
		#print(g.nodes.get(node_name).successors())
		for i in G.successors(node_name):
			print('\t'+str(i))
		print('predecessors:')
		for i in G.predecessors(node_name):
			print('\t'+str(i))


#### reactions

reactions_only_successors = []
for node_name in G.nodes.keys():
	node = G.nodes.get(node_name)
	if node['type']=='reaction':
		print('------------ '+str(node_name)+' ------------')
		print('successors:')
		#print(g.nodes.get(node_name).successors())
		for i in G.successors(node_name):
			print('\t'+str(i))
		print('predecessors:')
		for i in G.predecessors(node_name):
			print('\t'+str(i))


 #####

only_consumed_species = []
only_produced_species = []
for node_name in G.nodes.keys():
	node = G.nodes.get(node_name)
	if node['type']=='species':
		print('------------ '+str(node_name)+' ------------')
		if len(list(G.successors(node_name)))==0:
			if len(list(G.predecessors(node_name)))>0:
				only_produced_species.append(node_name)
				print('only produced')
		else:
			if len(list(G.predecessors(node_name)))==0:
				if len(list(G.successors(node_name)))>0:
					only_consumed_species.append(node_name)
					print('only_consumed')


print(only_consumed_species)
print(only_produced_species)


while


for species in only_consumed_species:
	#find the previous reactions of the species
	for node_name in G.predecessors(species):
		node = G.nodes.get(node_name)
		if node['type']=='reaction':
			#find if the reaction count is the size of the 
			for G.predecessors()




	node = G.nodes.get(species)






###
















rp_dict = rpsbml.readRPspecies(pathway_id)
reactions = list(rp_dict.keys())















g = nx.DiGraph()

#add the species
for species in model.getListOfSpecies():
	g.add_node(species.getId(),
			   type = 'species', 
			   label = species.getName())

#add the reactions
for reaction in model.getListOfReactions():
	g.add_node(reaction.getId(),
			   type = 'reaction',
			   label = reaction.getId())

for reaction in model.getListOfReactions():
	for reac in reaction.getListOfReactants():
		g.add_edge(reac.species, reaction.getId())
	for prod in reaction.getListOfProducts():
		g.add_edge(reaction.getId(), prod.species)


#### species


for node_name in g.nodes.keys():
	print('------------ '+str(node_name)+' ------------')
	node = g.nodes.get(node_name)
	if node['type']=='species':
		print('successors:')
		#print(g.nodes.get(node_name).successors())
		for i in g.successors(node_name):
			print('\t'+str(i))
		print('predecessors:')
		for i in g.predecessors(node_name):
			print('\t'+str(i))


#### reactions


for node_name in g.nodes.keys():
	print('------------ '+str(node_name)+' ------------')
	node = g.nodes.get(node_name)
	if node['type']=='reaction':
		print('successors:')
		#print(g.nodes.get(node_name).successors())
		for i in g.successors(node_name):
			print('\t'+str(i))
		print('predecessors:')
		for i in g.predecessors(node_name):
			print('\t'+str(i))


##### 


only_consumed_species = []
only_produced_species = []
for node_name in G.nodes.keys():
	print('------------ '+str(node_name)+' ------------')
	node = G.nodes.get(node_name)
	if node['type']=='species':
		if len(list(G.successors(node_name)))==0:
			if len(list(G.predecessors(node_name)))>0:
				only_produced_species.append(node_name)
				print('only produced')
		else:
			if len(list(G.predecessors(node_name)))==0:
				if len(list(G.successors(node_name)))>0:
					only_consumed_species.append(node_name)
					print('only_consumed')


print(only_consumed_species)
print(only_produced_species)









is_reac = True
reaction = None
reac_count = 1
for node_name in G.predecessors(start_species):
    node = G.nodes.get(node_name)
    if node['type']=='reaction':
        if reaction:
            print('Should only have one reaction predecessor for each species')
        else:
            reaction = node_name
    elif node['type']=='species':
        print('Species '+str(start_species)+' has species for predecessor')
    else:	
        print('Detected another node type: '+str(node['type']))



if reaction==None:
    print('Species '+str(start_species)+' does not have reaction predecessor')


print(reaction)
print('============================')
print(reac_count)


while is_reac:
    for spe_name in G.predecessors(reaction):
        spe_node = G.nodes.get(spe_name)
        if spe_node['type']=='species':
            for reac_name in G.predecessors(spe_name):
                print(reac_name)
                print(len(list(G.predecessors(reac_name))))
                print([i in only_consumed_species for i in G.predecessors(reac_name)])
                if all([i in only_consumed_species for i in G.predecessors(reac_name)]):
                	is_reac = False
                reac_node = G.nodes.get(reac_name)
                if reac_node['type']=='reaction':
                    if not reac_name==reaction:
                        reaction = reac_name
                        print('+++ '+str(reaction))
                        reac_count += 1
                elif reac_node['type']=='species':
                    print('Species '+str(start_species)+' has species for predecessor')
                else:
                    print('Detected another node type: '+str(node['type']))
        elif node['type']=='reaction':
            print('Reation '+str(spe_name)+' contains species as predecessor')
        else:
            print('Detected another node type: '+str(node['type']))


print(reac_count)
print('============================')

	