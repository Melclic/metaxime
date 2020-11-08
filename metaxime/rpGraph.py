import networkx as nx
from networkx.readwrite import json_graph
import logging
import os
import itertools
import numpy as np
import random

from .rpSBML import rpSBML


__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = []
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Melchior du Lac"
__status__ = "Development"


logging.basicConfig(
    level=logging.DEBUG,
    #level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

class rpGraph(rpSBML):
    """The class that hosts the networkx related functions
    """
    def __init__(self,
                 is_gem_sbml=False,
                 pathway_id='rp_pathway',
                 central_species_group_id='central_species',
                 sink_species_group_id='rp_sink_species',
                 model_name=None,
                 document=None,
                 path=None,
                 rpcache=None):
        """Constructor of the class

        Automatically constructs the network when calling the construtor

        :param rpcache: rpCache object
        :param model_name: The name of the model
        :param document: The libSBML Document of the model
        :param path: The path to the SBMKL file
        :param is_gem_sbml: Determine if its a full GEM model or not
        :param pathway_id: The Groups id of the heterologous pathway
        :param central_specues_group_id: The Groups id of the central species
        :param species_group_id: The id of the central species

        :type rpsbml: rpSBML
        :type pathway_id: str
        :type species_group_id: str

        .. document private functions
        .. automethod:: _makecomparegraphs
        .. automethod:: _makeGraph
        .. automethod:: _recursiveReacSuccessors
        .. automethod:: _recursiveReacPredecessors
        """
        super().__init__(model_name,
                         document,
                         path,
                         rpcache)
        self.logger = logging.getLogger(__name__)
        #WARNING: change this to reflect the different debugging levels
        self.logger.debug('Started instance of rpGraph')
        self.pathway_id = pathway_id
        self.central_species_group_id = central_species_group_id
        self.sink_species_group_id = sink_species_group_id
        self.G = None
        self.pathway_id = pathway_id
        self.num_reactions = 0
        self.num_species = 0
        if self.model:
            self._makeGraph(is_gem_sbml, pathway_id, central_species_group_id, sink_species_group_id)


    ######################################################################################################
    ######################################### Private Function ###########################################
    ######################################################################################################


    #TODO: if there are multiple EC number and multiple inchikeys, then you should construct all the alternative graphs and compare the, with your target, and return the one that has the highest score. Only then can you 
    def _makeCompareGraphs(self, inchikey_layers=2, ec_layers=3, pathway_id='rp_pathway'):
        """Make a special graphs for comparison whit the ID's being unique to the nodes

        Because comparisong of networkx graphs cannot use the attributes of the nodes, we create graphs based on the EC number of the reactions and the InChiKeys of the species

        :param inchikey_layers: The number of inchikey layers to take into consideration when comapring species
        :param ec_layers: The level of the EC number to take into consideration when comparing reactions
        :paran pathway_id: The heterologous pathway Groups id

        :rtype:
        :return:
        """
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


    #TODO: add the compartments to the species and reactions node descriptions
    def _makeGraph(self, is_gem_sbml=False, pathway_id='rp_pathway', central_species_group_id='central_species', sink_species_group_id='rp_sink_species'):
        """Private function that constructs the networkx graph

        :param is_gem_sbml: Determine what type of graph to build. If True then all the species and reactions will be added and not just the heterologous pathway.
        :param pathway_id: The pathway id of the heterologous pathway
        :param species_group_id: The id of the central species

        :type pathway_id: str
        :type species_group_id: str

        :return: None
        :rtype: None
        """
        #rp_species = [self.model.getSpecies(i) for i in self.readUniqueRPspecies(pathway_id)]
        groups = self.model.getPlugin('groups')
        self._checklibSBML(groups, 'Getting the groups')
        try:
            c_s = groups.getGroup(central_species_group_id)
            self._checklibSBML(c_s, 'Getting the annotation: '+str(central_species_group_id))
            rp_central_species_id = [i.getIdRef() for i in c_s.getListOfMembers()]
        except AttributeError:
            rp_central_species_id = []
        self.logger.debug('rp_central_species_id: '+str(rp_central_species_id))
        try:
            s_s = groups.getGroup(sink_species_group_id)
            self._checklibSBML(s_s, 'Getting the annotation: '+str(sink_species_group_id))
            rp_sink_species_id = [i.getIdRef() for i in s_s.getListOfMembers()]
        except AttributeError:
            rp_sink_species_id = []
        self.logger.debug('rp_sink_species_id: '+str(rp_sink_species_id))
        try:
            rp_pathway = groups.getGroup(pathway_id)
            rp_species_id = self.readUniqueRPspecies(pathway_id)
            rp_reactions_id = [i.getIdRef() for i in rp_pathway.getListOfMembers()]
            rp_annotation = self.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
        except AttributeError:
            rp_species_id = []
            rp_reactions_id = []
            rp_annotation = {}
        self.logger.debug('rp_species_id: '+str(rp_species_id))
        self.logger.debug('rp_reactions_id: '+str(rp_reactions_id))
        self.logger.debug('rp_annotation: '+str(rp_annotation))
        print('rp_reactions_id: '+str(rp_reactions_id))
        self.G = nx.DiGraph(brsynth=rp_annotation)
        #### add ALL the species and reactions ####
        #nodes
        for species in self.model.getListOfSpecies():
            is_central = False
            is_sink = False
            is_rp_pathway = False
            if species.getId() in rp_species_id:
                is_rp_pathway = True
            if species.getId() in rp_central_species_id:
                is_central = True
            if species.getId() in rp_sink_species_id:
                is_sink = True
            #add it if GEM then all, or if rp_pathway
            if is_rp_pathway or is_gem_sbml:
                self.num_species += 1
                self.G.add_node(species.getId(),
                                type='species',
                                name=species.getName(),
                                miriam=self.readMIRIAMAnnotation(species.getAnnotation()),
                                brsynth=self.readBRSYNTHAnnotation(species.getAnnotation()),
                                central_species=is_central,
                                sink_species=is_sink,
                                rp_pathway=is_rp_pathway)
        for reaction in self.model.getListOfReactions():
            is_rp_pathway = False
            if reaction.getId() in rp_reactions_id:
                is_rp_pathway = True
            if is_rp_pathway or is_gem_sbml:
                self.num_reactions += 1
            if reaction.getId() in rp_reactions_id or is_gem_sbml:
                self.G.add_node(reaction.getId(),
                                type='reaction',
                                miriam=self.readMIRIAMAnnotation(reaction.getAnnotation()),
                                brsynth=self.readBRSYNTHAnnotation(reaction.getAnnotation()),
                                rp_pathway=is_rp_pathway)
        #edges
        for reaction in self.model.getListOfReactions():
            if reaction.getId() in rp_reactions_id or is_gem_sbml:
                for reac in reaction.getListOfReactants():
                    self.G.add_edge(reac.species,
                                    reaction.getId(),
                                    stoichio=reac.stoichiometry)
                for prod in reaction.getListOfProducts():
                    self.G.add_edge(reaction.getId(),
                                    prod.species,
                                    stoichio=reac.stoichiometry)


    ## Recursive function that finds the order of the reactions in the graph 
    #
    # NOTE: only works for linear pathways... need to find a better way ie. Tree's
    #
    def _recursiveReacSuccessors(self, node_name, reac_list, all_res, num_reactions):
        """Recusrively order the reaction according to the successors

        :param node_name: The id of the starting node
        :param reac_list: The list of reactions that already have been run
        :param all_res:
        :param num_reactions:

        :type node_name: str
        :type reac_list: list
        :type all_res:
        :type num_reactions:

        :return: List of node ids
        :rtype: list
        """
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
        """Recusrively order the reaction according to the predecessors

        :param node_name: The id of the starting node
        :param reac_list: The list of reactions that already have been run
        :param all_res:
        :param num_reactions:

        :type node_name: str
        :type reac_list: list
        :type all_res:
        :type num_reactions:

        :return: List of node ids
        :rtype: list
        """
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



    '''
    def _recursiveHierarchy(self, node_name, num_nodes, ranked_nodes):
        self.G.successors(node_name)
    '''

    ######################################################################################################
    ########################################## Public Function ###########################################
    ######################################################################################################


    def onlyConsumedSpecies(self, only_central=False, only_rp_pathway=True):
        """Private function that returns the single parent species that are consumed only

        :param only_central: Focus on the central species only

        :type only_central: bool

        :return: List of node ids
        :rtype: list
        """
        only_consumed_species = []
        for node_name in self.G.nodes():
            node = self.G.node.get(node_name)
            if node['type']=='species':
                #NOTE: if central species then must also be rp_pathway species
                if (only_central and node['central_species']==True) or (only_rp_pathway and node['rp_pathway']==True) or (not only_central and not only_rp_pathway):
                    if not len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))==0:
                        only_consumed_species.append(node_name)
        return only_consumed_species


    def onlyProducedSpecies(self, only_central=False, only_rp_pathway=True):
        """Private function that returns the single parent produced species

        :param only_central: Focus on the central species only

        :type only_central: bool

        :return: List of node ids
        :rtype: list
        """
        only_produced_species = []
        for node_name in self.G.nodes():
            node = self.G.node.get(node_name)
            self.logger.debug('node_name: '+str(node_name))
            self.logger.debug('node: '+str(node))
            if node['type']=='species':
                #NOTE: if central species then must also be rp_pathway species
                if (only_central and node['central_species']==True) or (only_rp_pathway and node['rp_pathway']==True) or (not only_central and not only_rp_pathway):
                    if len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))>0:
                        only_produced_species.append(node_name)
        return only_produced_species


    ## Compare two rpgraph hypergraphs and return a score using a simple walk
    #
    # NOTE: source and target are used purely for clarity, you can inverse the two and the results are the same
    @staticmethod
    def compare(target_rpgraph, source_rpgraph, inchikey_layers=2, ec_layers=3, pathway_id='rp_pathway'):
        """Compare two rpgraph hypergraphs and return a score using a simple walk

        This function uses the gmatch4py package that implements a number of graph comparison 

        """
        import gmatch4py as gm
        #source_compare_graphs = source_rpgraph._makeCompareGraphs(inchikey_layers, ec_layers, pathway_id)
        source_compare_graphs = source_rpgraph._makeCompareGraphs(inchikey_layers, ec_layers, pathway_id)
        target_compare_graphs = target_rpgraph._makeCompareGraphs(inchikey_layers, ec_layers, pathway_id)
        #NOTE: here we use the greedy edit distance method but others may be used... 
        ged = gm.GreedyEditDistance(1,1,1,1)
        #ged = gm.GraphEditDistance(1,1,1,1)
        result = ged.compare([i[0] for i in source_compare_graphs]+[i[0] for i in target_compare_graphs], None)
        logging.debug('result: \n'+str([list(i) for i in result]))
        #seven is an arbitrary unit where 7==full info on ec number, smiles, inchikey etc...
        weights = np.array([sum(i[1])/7.0 for i in source_compare_graphs]+[sum(i[1])/7.0 for i in target_compare_graphs])
        weighted_similarity = np.array([i*weights for i in ged.similarity(result)])
        logging.debug('weighted_similarity: \n'+str([list(i) for i in weighted_similarity]))
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
        logging.debug('filtered_weighted_similarity: \n'+str([list(i) for i in filtered_weighted_similarity]))
        return max(map(max, filtered_weighted_similarity))


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
        for cons_cent_spe in self.onlyConsumedSpecies():
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
        for prod_cent_spe in self.onlyProducedSpecies():
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
