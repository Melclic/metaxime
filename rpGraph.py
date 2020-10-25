import networkx as nx
from networkx.readwrite import json_graph
import logging
import os
import itertools
import numpy as np
import random

from rpSBML import rpSBML

class rpGraph(rpSBML):
    """The class that hosts the networkx related functions
    """
    def __init__(self, pathway_id='rp_pathway', central_species_group_id='central_species', sink_species_group_id='rp_sink_species'):
        """Constructor of the class

        Automatically constructs the network when calling the construtor

        :param pathway_id: The pathway id of the heterologous pathway
        :param species_group_id: The id of the central species
        :param sink_species_group_id: The ids of the sink species

        :type pathway_id: str
        :type species_group_id: str
        :type sink_species_group_id: str

        .. document private functions
        .. automethod:: _makeCompareGraphs
        """
        super().__init__()
        self.logger = logging.getLogger(__name__)
        #WARNING: change this to reflect the different debugging levels
        self.logger.debug('Started instance of rpGraph')
        self.pathway_id = pathway_id
        self.central_species_group_id = central_species_group_id
        self.sink_species_group_id = sink_species_group_id
        self.G = None
        self.species = None
        self.reactions = None
        self.pathway_id = pathway_id
        self.num_reactions = 0
        self.central_species = []
        self.sink_species = []
        self.num_species = 0
        if rpsbml:
            self._makeGraph(pathway_id, central_species_group_id, sink_species_group_id)


    ######################################################################################################
    ######################################### Private Function ###########################################
    ######################################################################################################


    ## Make a special graphs for comparison whit the ID's being unique to the nodes
    #
    # Because comparisong of networkx graphs cannot use the attributes of the nodes, we create graphs based on the EC number 
    # of the reactions and the InChiKeys of the species
    #
    # TODO: if there are multiple EC number and multiple inchikeys, then you should construct all the alternative graphs and
    # compare the, with your target, and return the one that has the highest score. Only then can you 
    def _makeCompareGraphs(self, inchikey_layers=2, ec_layers=3, pathway_id='rp_pathway'):
        #retreive the pathway species and reactions
        species = [self.model.getSpecies(i) for i in self.readUniqueRPspecies(pathway_id)]
        groups = self.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        reactions = [self.model.getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]
        Gs = []
        #what you want to do is build a series of graphs that have the most info to the least 
        #WARNING: Probably need to checkt that there are no 2 species that have the same inchi_keys at a particular layer
        ###################### Species #######################
        # The species either have their inchikey's or their id's used as name of nodes. If smaller layers are input
        # then the full inchikey and their lower layers are input to build the graphs
        spe_comb = []
        spe_comb_info = []
        if inchikey_layers in [1,2,3]:
            for inch_lay in reversed(range(inchikey_layers, 4)):
                speid_newid = {}
                ############# InChiKey ###############
                for spe in species:
                    brsynth = self.readBRSYNTHAnnotation(spe.getAnnotation())
                    #loop through all the different layers of the inchikey and calculate the graph
                    if 'inchikey' in brsynth:
                        speid_newid[spe.getId()] = '-'.join(i for i in brsynth['inchikey'].split('-')[:inch_lay])
                    else:
                        miriam = self.readMIRIAMAnnotation(spe.getAnnotation())
                        if 'inchikey' in miriam:
                            if len(miriam['inchikey'])>1:
                                self.logger.warning('There are multiple inchikeys for '+str(spe.id)+': '+str(miriam['inchikey']))
                                self.logger.warning('Using the first one')
                            speid_newid[spe.getId()] = '-'.join(i for i in miriam['inchikey'][0].split('-')[:inch_lay])
                        else:
                            self.logger.warning('There is no inchikey associated with species: '+str(spe.getId()))
                            self.logger.warning('Setting species ID as the node id')
                            speid_newid[spe.getId()] = spe.getId()
                spe_comb.append(speid_newid)
                spe_comb_info.append(inch_lay)
        elif inchikey_layers==0:
            speid_newid = {}
            for spe in species:
                speid_newid[spe.getId()] = spe.getId()
            spe_comb.append(speid_newid)
            spe_comb_info.append(3) # ie. full info
        else:
            self.logger.error('Cannot recognise the inchi_layers input: '+str(inchikey_layers))
            return False
        ######################## Reactions ########################
        reac_comb = []
        reac_comb_info = []
        if ec_layers in [1,2,3,4]:
            for ec_lay in reversed(range(ec_layers, 5)):
                reacid_newid = {}
                ############### EC number ################
                for reac in reactions:
                    #brsynth = self.readBRSYNTHAnnotation(reac.getAnnotation())
                    miriam = self.readMIRIAMAnnotation(reac.getAnnotation())
                    self.logger.debug(str(reac)+' MIRIAM: '+str(miriam))
                    reacid_newid[reac.getId()] = []
                    if 'ec-code' in miriam:
                        #WARNING: need to deal with multiple ec numbers....
                        for ec_n in miriam['ec-code']:
                            reacid_newid[reac.getId()].append('.'.join(i for i in ec_n.split('.')[:ec_lay] if not i=='-')) #remove the layers that have '-' characters
                    else:
                        self.logger.warning('There is no EC number associated with reaction: '+str(reac.getId())+'. Setting the is as the node name')
                        reacid_newid[reac.getId()] = [reac.getId()] #consider random assignement of reaction id's since these may skew the comparison results
                reac_comb.append(reacid_newid)
                reac_comb_info.append(ec_lay)
        elif ec_layers==0:
            reacid_newid = {}
            reacid_newid[reac.getId()] = [reac.getId()]
            reac_comb.append(reacid_newid)
            reac_comb_info.append(4) #ie. full info
        elif ec_layers==-1:
            reacid_newid = {}
            reacid_newid[reac.getId()] = [brsynth['smiles'].upper()]
            reac_comb.append(reacid_newid)
            reac_comb_info.append(4) # ie. full info
        else:
            self.logger.error('Cannot interpret the ec_layers input: '+str(ec_layers))
            return False
        ###### make the graphs #####
        #remove the duplicates
        for rea in reac_comb:
            for rpx in rea:
                rea[rpx] = list(set(rea[rpx]))
        Gs = []
        #combine the different EC numbers per reaction
        #NOTE: These are ordered lists where the first values have the most info and the other have decreasing amounts
        self.logger.debug('spe_comb: '+str(spe_comb))
        self.logger.debug('spe_comb_info: '+str(spe_comb_info))
        self.logger.debug('reac_comb: '+str(reac_comb))
        self.logger.debug('reac_comb_info: '+str(reac_comb_info))
        self.logger.debug('----------------------------------------')
        for comb, comb_info in zip(list(itertools.product(spe_comb, reac_comb)), list(itertools.product(spe_comb_info, reac_comb_info))):
            self.logger.debug('comb: '+str(comb))
            self.logger.debug('comb_info: '+str(comb_info))
            ids = list(comb[1].keys())
            list_list = [comb[1][i] for i in ids]
            for r in list(itertools.product(*list_list)):
                tmp_reac = {key: value for (key, value) in zip(ids, r)}
                G = nx.DiGraph()
                for spe in species:
                    G.add_node(comb[0][spe.getId()])
                for rea in reactions:
                    G.add_node(tmp_reac[rea.getId()])
                for reaction in reactions:
                    for reac in reaction.getListOfReactants():
                        G.add_edge(comb[0][reac.species],
                                   tmp_reac[reaction.getId()])
                    for prod in reaction.getListOfProducts():
                        G.add_edge(tmp_reac[reaction.getId()],
                                   comb[0][prod.species])
                Gs.append((G, comb_info))
        return Gs


    ################################# Analyse and make graph #####################

    def _makeGraph(self, pathway_id='rp_pathway', central_species_group_id='central_species', sink_species_group_id='rp_sink_species'):
        """Private function that constructs the networkx graph

        :param pathway_id: The pathway id of the heterologous pathway
        :param species_group_id: The id of the central species

        :type pathway_id: str
        :type species_group_id: str

        :return: None
        :rtype: None
        """
        self.species = [self.model.getSpecies(i) for i in self.readUniqueRPspecies(pathway_id)]
        groups = self.model.getPlugin('groups')
        c_s = groups.getGroup(central_species_group_id)
        self.central_species = [i.getIdRef() for i in c_s.getListOfMembers()]
        s_s = groups.getGroup(sink_species_group_id)
        self.sink_species = [i.getIdRef() for i in s_s.getListOfMembers()]
        rp_pathway = groups.getGroup(pathway_id)
        self.reactions = [self.model.getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]
        self.G = nx.DiGraph(brsynth=self.readBRSYNTHAnnotation(rp_pathway.getAnnotation()))
        #nodes
        for spe in self.species:
            self.num_species += 1
            is_central = False
            is_sink = False
            if spe.getId() in self.central_species:
                is_central = True
            if spe.getId() in self.sink_species:
                is_sink = True
            self.G.add_node(spe.getId(),
                            type='species',
                            name=spe.getName(),
                            miriam=self.readMIRIAMAnnotation(spe.getAnnotation()),
                            brsynth=self.readBRSYNTHAnnotation(spe.getAnnotation()),
                            central_species=is_central,
                            sink_species=is_sink)
        for reac in self.reactions:
            self.num_reactions += 1
            self.G.add_node(reac.getId(),
                            type='reaction',
                            miriam=self.readMIRIAMAnnotation(reac.getAnnotation()),
                            brsynth=self.readBRSYNTHAnnotation(reac.getAnnotation()))
        #edges
        for reaction in self.reactions:
            for reac in reaction.getListOfReactants():
                self.G.add_edge(reac.species,
                                reaction.getId(),
                                stoichio=reac.stoichiometry)
            for prod in reaction.getListOfProducts():
                self.G.add_edge(reaction.getId(),
                                prod.species,
                                stoichio=reac.stoichiometry)


    def _onlyConsumedSpecies(self, only_central=False):
        """Private function that returns the single parent species that are consumed only

        :return: List of node ids
        :rtype: list
        """
        only_consumed_species = []
        for node_name in self.G.nodes():
            node = self.G.node.get(node_name)
            if node['type']=='species':
                if not len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))==0:
                    only_consumed_species.append(node_name)
        return only_consumed_species


    def _onlyProducedSpecies(self, only_central=False):
        """Private function that returns the single parent produced species

        :return: List of node ids
        :rtype: list
        """
        only_produced_species = []
        for node_name in self.G.nodes():
            node = self.G.node.get(node_name)
            self.logger.debug('node_name: '+str(node_name))
            self.logger.debug('node: '+str(node))
            if node['type']=='species':
                if only_central:
                    if node['central_species']==True:
                        if len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))>0:
                            only_produced_species.append(node_name)
                else:
                    if len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))>0:
                        only_produced_species.append(node_name)
        return only_produced_species


    ## Recursive function that finds the order of the reactions in the graph 
    #
    # NOTE: only works for linear pathways... need to find a better way ie. Tree's
    #
    def _recursiveReacSuccessors(self, node_name, reac_list, all_res, num_reactions):
        current_reac_list = [i for i in reac_list]
        self.logger.debug('-------- '+str(node_name)+' --> '+str(reac_list)+' ----------')
        succ_node_list = [i for i in self.G.successors(node_name)]
        flat_reac_list = [i for sublist in reac_list for i in sublist]
        self.logger.debug('flat_reac_list: '+str(flat_reac_list))
        self.logger.debug('current_reac_list: '+str(current_reac_list))
        if len(flat_reac_list)==num_reactions:
            self.logger.debug('Returning')
            #return current_reac_list
            all_res.append(current_reac_list)
            return all_res
        self.logger.debug('succ_node_list: '+str(succ_node_list))
        if not succ_node_list==[]:
            #can be multiple reactions at a given step
            multi_reac = []
            for n_n in succ_node_list:
                n = self.G.node.get(n_n)
                if n['type']=='reaction':
                    if not n_n in flat_reac_list:
                        multi_reac.append(n_n)
            #remove the ones that already exist
            self.logger.debug('multi_reac: '+str(multi_reac))
            multi_reac = [x for x in multi_reac if x not in flat_reac_list]
            self.logger.debug('multi_reac: '+str(multi_reac))
            if multi_reac:
                current_reac_list.append(multi_reac)
            self.logger.debug('current_reac_list: '+str(current_reac_list))
            #loop through all the possibilities
            for n_n in succ_node_list:
                n = self.G.node.get(n_n)
                if n['type']=='reaction':
                    if n_n in multi_reac:
                        self._recursiveReacSuccessors(n_n, current_reac_list, all_res, num_reactions)
                elif n['type']=='species':
                    if n['central_species']==True:
                        self._recursiveReacSuccessors(n_n, current_reac_list, all_res, num_reactions)
        return all_res


    ##
    #
    # NOTE: only works for linear pathways... need to find a better way
    #
    def _recursiveReacPredecessors(self, node_name, reac_list, all_res, num_reactions):
        current_reac_list = [i for i in reac_list]
        self.logger.debug('-------- '+str(node_name)+' --> '+str(reac_list)+' ----------')
        pred_node_list = [i for i in self.G.predecessors(node_name)]
        flat_reac_list = [i for sublist in reac_list for i in sublist]
        self.logger.debug('flat_reac_list: '+str(flat_reac_list))
        self.logger.debug('current_reac_list: '+str(current_reac_list))
        if len(flat_reac_list)==num_reactions:
            self.logger.debug('Returning')
            #return current_reac_list
            all_res.append(current_reac_list)
            return all_res

            #can be multiple reactions at a given step
            multi_reac = []
            for n_n in pred_node_list:
                n = self.G.node.get(n_n)
                if n['type']=='reaction':
                    if not n_n in flat_reac_list:
                        multi_reac.append(n_n)
            #remove the ones that already exist
            self.logger.debug('multi_reac: '+str(multi_reac))
            multi_reac = [x for x in multi_reac if x not in flat_reac_list]
            self.logger.debug('multi_reac: '+str(multi_reac))
            if multi_reac:
                current_reac_list.append(multi_reac)
            self.logger.debug('current_reac_list: '+str(current_reac_list))
            #loop through all the possibilities
            for n_n in pred_node_list:
                n = self.G.node.get(n_n)
                if n['type']=='reaction':
                    if n_n in multi_reac:
                        self._recursiveReacPredecessors(n_n, current_reac_list, all_res, num_reactions)
                elif n['type']=='species':
                    if n['central_species']==True:
                        self._recursiveReacPredecessors(n_n, current_reac_list, all_res, num_reactions)
        return all_res


    '''
    def _recursiveHierarchy(self, node_name, num_nodes, ranked_nodes):
        self.G.successors(node_name)
    '''

    ######################################################################################################
    ########################################## Public Function ###########################################
    ######################################################################################################


    ## Compare two rpgraph hypergraphs and return a score using a simple walk
    #
    # NOTE: source and target are used purely for clarity, you can inverse the two and the results are the same
    def compare(source_rpgraph, target_rpgraph, inchikey_layers=2, ec_layers=3, pathway_id='rp_pathway'):
        import gmatch4py as gm
        source_compare_graphs = source_rpgraph._makeCompareGraphs(inchikey_layers, ec_layers, pathway_id)
        target_compare_graphs = target_rpgraph._makeCompareGraphs(inchikey_layers, ec_layers, pathway_id)
        #NOTE: here we use the greedy edit distance method but others may be used... 
        ged = gm.GreedyEditDistance(1,1,1,1)
        #ged = gm.GraphEditDistance(1,1,1,1)
        result = ged.compare([i[0] for i in source_compare_graphs]+[i[0] for i in target_compare_graphs], None)
        rpGraph.logger.debug('result: \n'+str([list(i) for i in result]))
        #seven is an arbitrary unit where 7==full info on ec number, smiles, inchikey etc...
        weights = np.array([sum(i[1])/7.0 for i in source_compare_graphs]+[sum(i[1])/7.0 for i in target_compare_graphs])
        weighted_similarity = np.array([i*weights for i in ged.similarity(result)])
        rpGraph.logger.debug('weighted_similarity: \n'+str([list(i) for i in weighted_similarity]))
        #weighted_distance = np.array([i*weights for i in ged.distance(result)])
        #rpGraph.logger.debug('weighted_distance: \n'+str([list(i) for i in weighted_distance]))
        filtered_weighted_similarity =  []
        source_pos = [i for i in range(len(source_compare_graphs))]
        for i in range(len(weighted_similarity)):
            tmp = []
            for y in range(len(weighted_similarity[i])):
                if i in source_pos and not y in source_pos:
                    tmp.append(weighted_similarity[i][y])
                elif i not in source_pos and y in source_pos:
                    tmp.append(weighted_similarity[i][y])
                else:
                    tmp.append(0.0)
            filtered_weighted_similarity.append(tmp)
        rpGraph.logger.debug('filtered_weighted_similarity: \n'+str([list(i) for i in filtered_weighted_similarity]))
        return max(map(max, filtered_weighted_similarity))


    def _recursiveReacPredecessors(self, node_name, reac_list):
        """Return the next linear predecessors

        Better than before, however bases itself on the fact that central species do not have multiple predesessors
        if that is the case then the algorithm will return badly ordered reactions

        :param node_name: The id of the starting node
        :param reac_list: The list of reactions that already have been run

        :type node_name: str
        :type reac_list: list

        :return: List of node ids
        :rtype: list
        """
        self.logger.debug('-------- '+str(node_name)+' --> '+str(reac_list)+' ----------')
        pred_node_list = [i for i in self.G.predecessors(node_name)]
        self.logger.debug(pred_node_list)
        if pred_node_list==[]:
            return reac_list
        for n_n in pred_node_list:
            n = self.G.node.get(n_n)
            if n['type']=='reaction':
                if n_n in reac_list:
                    return reac_list
                else:
                    reac_list.append(n_n)
                    self._recursiveReacPredecessors(n_n, reac_list)
            elif n['type']=='species':
                if n['central_species']==True:
                    self._recursiveReacPredecessors(n_n, reac_list)
                else:
                    return reac_list
        return reac_list


    def orderedRetroReactions(self):
        """Public function to return the linear list of reactions

        :return: List of node ids
        :rtype: list
        """
        #Note: may be better to loop tho
        for prod_spe in self._onlyProducedSpecies():
            self.logger.debug('Testing '+str(prod_spe))
            ordered = self._recursiveReacPredecessors(prod_spe, [])
            self.logger.debug(ordered)
            if len(ordered)==self.num_reactions:
                return [i for i in reversed(ordered)]
        self.logger.error('Could not find the full ordered reactions')
        return []

    ############################# graph analysis ################################


    def exportJSON(self):
        return json_graph.node_link_data(self.G)


    ## Warning that this search algorithm only works for mono-component that are not networks (i.e where reactions follow each other)
    # DEPRECATED: this is linear
    # NOTE: only works for linear pathways... need to find a better way
    #
    def orderedRetroReactions(self):
        #Note: may be better to loop tho
        succ_res = []
        for cons_cent_spe in self._onlyConsumedSpecies():
            res = self._recursiveReacSuccessors(cons_cent_spe, [], [], self.num_reactions)
            if res:
                self.logger.debug(res)
                if len(res)==1:
                    succ_res = res[0]
                else:
                    self.logger.error('Multiple successors results: '+str(res))
            else:
                self.logger.warning('Successors no results')
        prod_res = []
        for prod_cent_spe in self._onlyProducedSpecies():
            res = self._recursiveReacPredecessors(prod_cent_spe, [], [], self.num_reactions)
            if res:
                self.logger.debug(res)
                if len(res)==1:
                    prod_res = [i for i in reversed(res[0])]
                else:
                    self.logger.error('Mutliple predecessors results: '+str(res))
            else:
                self.logger.warning('Predecessors no results')
        if succ_res and prod_res:
            if not succ_res==prod_res:
                self.logger.warning('Both produce results and are not the same')
                self.logger.warning('succ_res: '+str(succ_res))
                self.logger.warning('prod_res: '+str(prod_res))
            else:
                self.logger.debug('Found solution: '+str(succ_res))
                return succ_res
        return []


