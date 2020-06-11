import networkx as nx
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


##
#
#
class rpGraph:
    ##
    #
    #
    def __init__(self, rpsbml, pathway_id='rp_pathway', species_group_id='central_species'):
        self.rpsbml = rpsbml
        self.logger = logging.getLogger(__name__)
        #WARNING: change this to reflect the different debugging levels
        self.logger.info('Started instance of rpGraph')
        self.pathway_id = pathway_id
        self.species_group_id = species_group_id
        self.G = None
        self.species = None
        self.reactions = None
        self.pathway_id = pathway_id
        self.num_reactions = 0
        self.num_species = 0
        self._makeGraph(pathway_id, species_group_id)


    ##
    #
    #
    def _makeGraph(self, pathway_id='rp_pathway', species_group_id='central_species'):
        self.species = [self.rpsbml.model.getSpecies(i) for i in self.rpsbml.readUniqueRPspecies(pathway_id)]
        groups = self.rpsbml.model.getPlugin('groups')
        central_species = groups.getGroup(species_group_id)
        central_species = [i.getIdRef() for i in central_species.getListOfMembers()]
        rp_pathway = groups.getGroup(pathway_id)
        self.reactions = [self.rpsbml.model.getReaction(i.getIdRef()) for i in rp_pathway.getListOfMembers()]
        self.G = nx.DiGraph(brsynth=self.rpsbml.readBRSYNTHAnnotation(rp_pathway.getAnnotation()))
        #nodes
        for spe in self.species:
            self.num_species += 1
            is_central = False
            if spe.getId() in central_species:
                is_central = True
            self.G.add_node(spe.getId(), 
                            type='species',
                            name=spe.getName(),
                            miriam=self.rpsbml.readMIRIAMAnnotation(spe.getAnnotation()),
                            brsynth=self.rpsbml.readBRSYNTHAnnotation(spe.getAnnotation()),
                            central_species=is_central)
        for reac in self.reactions:
            self.num_reactions += 1
            self.G.add_node(reac.getId(),
                            type='reaction',
                            miriam=self.rpsbml.readMIRIAMAnnotation(reac.getAnnotation()),
                            brsynth=self.rpsbml.readBRSYNTHAnnotation(reac.getAnnotation()))
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


    ##
    #
    #
    def _onlyConsumedSpecies(self):
        only_consumed_species = []
        for node_name in self.G.nodes.keys():
            node = self.G.nodes.get(node_name)
            if node['type']=='species':
                if not len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))==0:
                    only_consumed_species.append(node_name)
        return only_consumed_species


    ##
    #
    #
    def _onlyConsumedCentralSpecies(self):
        only_consumed_species = []
        for node_name in self.G.nodes.keys():
            node = self.G.nodes.get(node_name)
            if node['type']=='species':
                if node['central_species']==True:
                    if not len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))==0:
                        only_consumed_species.append(node_name)
        return only_consumed_species



    ##
    #
    #
    def _onlyProducedSpecies(self):
        only_produced_species = []
        for node_name in self.G.nodes.keys():
            node = self.G.nodes.get(node_name)
            if node['type']=='species':
                if len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))>0:
                    only_produced_species.append(node_name)
        return only_produced_species


    ##
    #
    #
    def _onlyProducedCentralSpecies(self):
        only_produced_species = []
        for node_name in self.G.nodes.keys():
            node = self.G.nodes.get(node_name)
            if node['type']=='species':
                if node['central_species']==True:
                    if len(list(self.G.successors(node_name)))==0 and len(list(self.G.predecessors(node_name)))>0:
                        only_produced_species.append(node_name)
        return only_produced_species


    ## 
    #
    # Better than before, however bases itself on the fact that central species do not have multiple predesessors
    # if that is the case then the algorithm will return badly ordered reactions
    #
    def _recursiveReacPrecessors(self, node_name, reac_list):
        self.logger.info('-------- '+str(node_name)+' --> '+str(reac_list)+' ----------')
        pred_node_list = [i for i in self.G.predecessors(node_name)]
        self.logger.info(pred_node_list)
        if pred_node_list==[]:
            return reac_list
        for n_n in pred_node_list:
            n = self.G.nodes.get(n_n)
            if n['type']=='reaction':
                if n_n in reac_list:
                    return reac_list
                else:
                    reac_list.append(n_n)
                    self._recursiveReacPrecessors(n_n, reac_list)
            elif n['type']=='species':
                if n['central_species']==True:
                    self._recursiveReacPrecessors(n_n, reac_list)
                else:
                    return reac_list
        return reac_list


    ## Warning that this search algorithm only works for mono-component reactions
    #
    #
    def orderedRetroReaction(self):
        #Note: may be better to loop tho
        if not len(prod_cent)==1:
            self.logger.warning('The rppathway contains more than one only produced central species: '+str(prod_cent))
        for i in prod_cent:
            self.logger.info('Testing '+str(i))
            ordered = self._recursiveReacPrecessors(i, [])
            self.logger.info(ordered)
            if len(ordered)==self.num_reactions:
                return ordered
            

    ################################################# BELOW IS DEV ################################

    """
    def orderedRetroReac(self):
        for node_name in self.G.nodes:
            node = self.G.nodes.get(node_name)
            if node['type']=='reaction':
                self.logger.info('---> Starting reaction: '+str(node_name))
                tmp_retrolist = [node_name]
                is_not_end_reac = True
                while is_not_end_reac:
                    tmp_spereacs = []
                    #for all the species, gather the ones that are not valid
                    for spe_name in self.G.predecessors(tmp_retrolist[-1]):
                        spe_node = self.G.nodes.get(spe_name)
                        if spe_node['type']=='reaction':
                            self.logger.warning('Reaction '+str(tmp_retrolist[-1])+' is directly connected to the following reaction '+str(spe_name))
                            continue
                        elif spe_node['type']=='species':
                            if spe_node['central_species']==True:
                                self.logger.info('\tSpecies: '+str(spe_name))
                                self.logger.info('\t'+str([i for i in self.G.predecessors(spe_name)]))
                                tmp_spereacs.append([i for i in self.G.predecessors(spe_name)])
                        else:
                            self.logger.warning('Node type should be either reaction or species: '+str(node['type']))
                    #remove empty lists
                    self.logger.info(tmp_spereacs)
                    tmp_spereacs = [i for i in tmp_spereacs if i != []]
                    self.logger.info(tmp_spereacs)
                    #return the number of same intersect
                    if tmp_spereacs==[]:
                        is_not_end_reac = False
                        continue
                    tmp_spereacs = list(set.intersection(*map(set, tmp_spereacs)))
                    self.logger.info(tmp_spereacs)
                    if len(tmp_spereacs)>1:     
                        self.logger.warning('There are multiple matches: '+str(tmp_spereacs))
                    elif len(tmp_spereacs)==0:
                        self.logger.info('Found the last reaction')
                        is_not_end_reac = False
                    elif len(tmp_spereacs)==1:
                        self.logger.info('Found the next reaction: '+str(tmp_spereacs[0]))
                        if tmp_spereacs[0] not in tmp_retrolist:
                            tmp_retrolist.append(tmp_spereacs[0])
                        else:
                            self.logger.warning('Trying to add a reaction in the sequence that already exists')
                            is_not_end_reac = False
                    self.logger.info(tmp_retrolist)
                self.logger.info('The tmp result is: '+str(tmp_retrolist))
                if len(tmp_retrolist)==self.num_reactions:
                    return tmp_retrolist


    def orderedReac(self):
        for node_name in self.G.nodes:
            node = self.G.nodes.get(node_name)
            if node['type']=='reaction':
                self.logger.info('---> Starting reaction: '+str(node_name))
                tmp_retrolist = [node_name]
                is_not_end_reac = True
                while is_not_end_reac:
                    tmp_spereacs = []
                    #for all the species, gather the ones that are not valid
                    for spe_name in self.G.successors(tmp_retrolist[-1]):
                        spe_node = self.G.nodes.get(spe_name)
                        if spe_node['type']=='reaction':
                            self.logger.warning('Reaction '+str(tmp_retrolist[-1])+' is directly connected to the following reaction '+str(spe_name))
                            continue
                        elif spe_node['type']=='species':
                            if spe_node['central_species']==True:
                                self.logger.info('\tSpecies: '+str(spe_name))
                                self.logger.info('\t'+str([i for i in self.G.successors(spe_name)]))
                                tmp_spereacs.append([i for i in self.G.successors(spe_name)])
                        else:
                            self.logger.warning('Node type should be either reaction or species: '+str(node['type']))
                    #remove empty lists
                    self.logger.info(tmp_spereacs)
                    tmp_spereacs = [i for i in tmp_spereacs if i!=[]]
                    self.logger.info(tmp_spereacs)
                    #return the number of same intersect
                    if tmp_spereacs==[]:
                        is_not_end_reac = False
                        continue
                    tmp_spereacs = list(set.intersection(*map(set, tmp_spereacs)))
                    self.logger.info(tmp_spereacs)
                    if len(tmp_spereacs)>1:     
                        self.logger.warning('There are multiple matches: '+str(tmp_spereacs))
                    elif len(tmp_spereacs)==0:
                        self.logger.info('Found the last reaction')
                        is_not_end_reac = False
                    elif len(tmp_spereacs)==1:
                        self.logger.info('Found the next reaction: '+str(tmp_spereacs[0]))
                        if tmp_spereacs[0] not in tmp_retrolist:
                            tmp_retrolist.append(tmp_spereacs[0])
                        else:
                            self.logger.warning('Trying to add a reaction in the sequence that already exists')
                            is_not_end_reac = False
                    self.logger.info(tmp_retrolist)
                self.logger.info('The tmp result is: '+str(tmp_retrolist))
                if len(tmp_retrolist)==self.num_reactions:
                    return tmp_retrolist

    """
