
species = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8']
reactions = ['R1', 'R2', 'R3', 'R4']

stoichio = [
	[1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	[0.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 0.0],
	[0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, -1.0],
	[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0]
]

stoichio_trans = [
	[1.0, 0.0, 0.0, 0.0],
	[-1.0, 1.0, 0.0, 0.0],
	[-1.0, 0.0, 0.0, 0.0],
	[0.0, -1.0, 1.0, 0.0],
	[0.0, 1.0, 0.0, 0.0],
	[0.0, 0.0, -1.0, 1.0],
	[0.0, 0.0, 0.0, -1.0],
	[0.0, 0.0, -1.0, 1.0]
]


## Function to find out the order of the reactions that connect a start and end node
#
#
def workupGraph(start_reaction, end_reaction):
	is_not_end = True
	next_reaction
	while is_not_end:
		count = 0
		for reac in stoichio:
			if start_reaction==reactions[count]:

			count += 1



count = 0
start_reaction = []
end_reaction = []
order_reactions = []
for reac in stoichio:
	print('=========== '+str(reactions[count])+' =====================')
	#################### consumed #################
	consumed_species = [species[i] for i in range(len(reac)) if reac[i]<0.0]
	print('consumed_species')
	print(consumed_species)
	print('produced_species')
	print(produced_species)
	consumed_only_species = []
	for spe in consumed_species:
		#check if consumed species is not produced anywhere else
		if stoichio_trans[species.index(spe)].count(1.0)==0:
			consumed_only_species.append(spe)
	print('consumed_only_species')
	print(consumed_only_species)
	if consumed_species==consumed_only_species and consumed_species and consumed_only_species:
		print('Consumed species is start')
		start_reaction.append(reactions[count])
	##################### produced ####################
	produced_species = [species[i] for i in range(len(reac)) if reac[i]>0.0]
	produced_only_species = []
	for spe in produced_species:
		#check that all the consumed species are not further connected
		if stoichio_trans[species.index(spe)].count(-1.0)==0:
			produced_only_species.append(spe)
	print('produced_only_species')
	print(produced_only_species)
	if produced_species==produced_only_species and produced_species and produced_only_species:
		print('Produced species is end')
		end_reaction.append(reactions[count])
	count += 1
	
print('start_reaction')
print(start_reaction)
print('end_reaction')
print(end_reaction)





for spe in stoichio_trans:
	#means either single produced or single consumed
	if sum([1.0 for x in column if x!=0.0])==1.0: