import rpSBML
import numpy as np

pathway_id = 'rp_pathway'

rpsbml = rpSBML.rpSBML('rp_4_5')
rpsbml.readSBML('rpglobalscore/rp_4_5.sbml.xml')

species = rpsbml.readUniqueRPspecies(pathway_id)
rp_dict = rpsbml.readRPspecies(pathway_id)
reactions = list(rp_dict.keys())


#stochio_matrix = [[0.0]*len(species)]*len(reactions)
matrix = []
stoichio = []
for i in range(len(reactions)):
	a = []
	b = []
	for y in range(len(species)):
		a.append(0.0)
		b.append(0.0)
	matrix.append(a)
	stoichio.append(b)


for i in matrix:
	print(i)	

for i in stoichio:
	print(i)	


for reac in rp_dict:
	#print(str(reac)+'-> '+str(reactions.index(reac)))
	#print('\treactants')
	for spe in rp_dict[reac]['reactants']:
		#print('\t\t'+str(spe)+' -> '+str(species.index(spe)))
		#print('\t\t['+str(reactions.index(reac))+', '+str(species.index(spe))+']')
		if not rp_dict[reac]['reactants'][spe]==0.0:
			matrix[reactions.index(reac)][species.index(spe)] = -1.0
		stoichio[reactions.index(reac)][species.index(spe)] = -rp_dict[reac]['reactants'][spe]
	#print('\tproducts:')
	for spe in rp_dict[reac]['products']:
		#print('\t\t'+str(spe)+' -> '+str(species.index(spe)))
		#print('\t\t['+str(reactions.index(reac))+', '+str(species.index(spe))+']')
		if not rp_dict[reac]['products'][spe]==0.0:
			matrix[reactions.index(reac)][species.index(spe)] = 1.0
		stoichio[reactions.index(reac)][species.index(spe)] = rp_dict[reac]['products'][spe]


for i in matrix:
	print(i)

for i in stoichio:
	print(i)


####################################################################################################
####################################################################################################
####################################################################################################


matrix_trans = []
stoichio_trans = []
for i in range(len(species)):
	a = []
	b = []
	for y in range(len(reactions)):
		a.append(0.0)
		b.append(0.0)
	matrix_trans.append(a)
	stoichio_trans.append(b)


for i in matrix_trans:
	print(i)


for i in stoichio_trans:
	print(i)


for reac in rp_dict:
	#print(str(reac)+'-> '+str(reactions.index(reac)))
	#print('\treactants')
	for spe in rp_dict[reac]['reactants']:
		#print('\t\t'+str(spe)+' -> '+str(species.index(spe)))
		#print('\t\t['+str(reactions.index(reac))+', '+str(species.index(spe))+']')
		if not rp_dict[reac]['reactants'][spe]==0.0:
			matrix_trans[species.index(spe)][reactions.index(reac)] = -1.0
		stoichio_trans[species.index(spe)][reactions.index(reac)] = -rp_dict[reac]['reactants'][spe]
	#print('\tproducts:')
	for spe in rp_dict[reac]['products']:
		#print('\t\t'+str(spe)+' -> '+str(species.index(spe)))
		#print('\t\t['+str(reactions.index(reac))+', '+str(species.index(spe))+']')
		if not rp_dict[reac]['products'][spe]==0.0:
			matrix_trans[species.index(spe)][reactions.index(reac)] = 1.0
		stoichio_trans[species.index(spe)][reactions.index(reac)] = rp_dict[reac]['products'][spe]


for i in matrix_trans:
	print(i)


for i in stoichio_trans:
	print(i)






#find the nodes at the end of the directed graph
end_species = []
for spe in species:
	column = []
	for reac in range(len(reactions)):
		column.append(stochio_matrix[reac][species.index(spe)])
	print(column)
	if sum([1.0 for x in column if x!=0.0])==1.0:
		try:
			print(stochio_matrix[column.index(1.0)])
			if sum(stochio_matrix[column.index(1.0)])==0.0:
				end_species.append(spe)
				print(spe)
		except ValueError:
			pass





#find the nodes at the start of the directed graph
start_species = []
for spe in species:
	column = []
	for reac in range(len(reactions)):
		column.append(stochio_matrix[reac][species.index(spe)])
	if sum([1.0 for x in column if x!=0.0])==1.0:
		try:
			print(stochio_matrix[column.index(-1.0)])
			if sum(stochio_matrix[column.index(-1.0)])==1.0:
				start_species.append(spe)
		except ValueError:
			pass






for reac in reactions:
	print(reac)
	reactions.index(reac)
	for spe in species:
		column = 
		if 











df2 = pd.DataFrame(np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]),
                   columns=['a', 'b', 'c'])











for row in stochio_matrix:








df2 = pd.DataFrame(stochio_matrix,
                   columns=species,
                   index=reactions)


















