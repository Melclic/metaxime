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
    def __init__(self, rpsbml, pathway_id='rp_pathway'):
        self.rpsbml = rpsbml
        self.logger = logging.getLogger(__name__)
        #WARNING: change this to reflect the different debugging levels
        self.logger.info('Started instance of rpSBML')
        self.pathway_id = pathway_id
        self.G = None
        self.species = None
        self.reactions = None
        self.pathway_id = pathway_id
        self.num_reactions = 0
        self.num_species = 0
        self._makeGraph(pathway_id)


    ##
    #
    #
    def _makeGraph(self, pathway_id='rp_pathway'):
        self.G = nx.DiGraph()
        self.species = [self.rpsbml.model.getSpecies(i) for i in self.rpsbml.readUniqueRPspecies(pathway_id)]
        rp_dict = self.rpsbml.readRPspecies(pathway_id)
        self.reactions = [self.rpsbml.model.getReaction(i) for i in list(rp_dict.keys())]
        for spe in self.species:
            self.num_species += 1
            self.G.add_node(spe.getId(), type='species')
        for reac in self.reactions:
            self.num_reactions += 1
            self.G.add_node(reac.getId(), type='reaction')
        for reaction in self.reactions:
            for reac in reaction.getListOfReactants():
                self.G.add_edge(reac.species, reaction.getId())
            for prod in reaction.getListOfProducts():
                self.G.add_edge(reaction.getId(), prod.species)


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
    def _countBackReac(self, start_species):
        is_reac = True
        reaction = None
        ordered_reac = []
        self.logger.info('============= '+str(start_species)+' ==============')
        only_consumed_species = self._onlyConsumedSpecies()
        for node_name in self.G.predecessors(start_species):
            self.logger.info('-> '+str(node_name))
            node = self.G.nodes.get(node_name)
            if node['type']=='reaction':
                if reaction:
                    self.logger.warning('Species is produced by multiple reactions')
                    return []
                else:
                    reaction = node_name
                    ordered_reac.append(node_name)
            elif node['type']=='species':
                self.logger.warning('Species '+str(start_species)+' has species for predecessor')
            else:
                self.logger.warning('Detected another node type: '+str(node['type']))
        if reaction==None:
            self.logger.error('Species '+str(start_species)+' does not have reaction predecessor')
            return 0
        while is_reac:
            for spe_name in self.G.predecessors(reaction):
                spe_node = self.G.nodes.get(spe_name)
                if spe_node['type']=='species':
                    for reac_name in self.G.predecessors(spe_name):
                        if all([i in only_consumed_species for i in self.G.predecessors(reac_name)]):
                            is_reac = False
                        reac_node = self.G.nodes.get(reac_name)
                        if reac_node['type']=='reaction':
                            if not reac_name==reaction:
                                reaction = reac_name
                                ordered_reac.append(reac_name)
                        elif reac_node['type']=='species':
                            self.logger.warning('Species '+str(start_species)+' has species for predecessor')
                        else:
                            self.logger.warning('Detected another node type: '+str(node['type']))
                elif node['type']=='reaction':
                    self.logger.warning('Reation '+str(spe_name)+' contains species as predecessor')
                else:
                    self.logger.warning('Detected another node type: '+str(node['type']))
        return ordered_reac


    ##
    #
    #
    def _countForwReac(self, end_species):
        is_reac = True
        reaction = None
        ordered_reac = []
        self.logger.info('============= '+str(end_species)+' ==============')
        only_produced_species = self._onlyProducedSpecies()
        for node_name in self.G.successors(end_species):
            self.logger.info('-> '+str(node_name))
            node = self.G.nodes.get(node_name)
            if node['type']=='reaction':
                if reaction:
                    self.logger.warning('Species is produced by multiple reactions')
                    return []
                else:
                    reaction = node_name
                    ordered_reac.append(node_name)
            elif node['type']=='species':
                self.logger.warning('Species '+str(end_species)+' has species for predecessor')
            else:
                self.logger.warning('Detected another node type: '+str(node['type']))
        if reaction==None:
            self.logger.error('Species '+str(end_species)+' does not have reaction predecessor')
            return 0
        while is_reac:
            for spe_name in self.G.successors(reaction):
                spe_node = self.G.nodes.get(spe_name)
                if spe_node['type']=='species':
                    for reac_name in self.G.successors(spe_name):
                        if all([i in only_produced_species for i in self.G.successors(reac_name)]):
                            is_reac = False
                        reac_node = self.G.nodes.get(reac_name)
                        if reac_node['type']=='reaction':
                            if not reac_name==reaction:
                                reaction = reac_name
                                ordered_reac.append(reac_name)
                        elif reac_node['type']=='species':
                            self.logger.warning('Species '+str(end_species)+' has species for predecessor')
                        else:
                            self.logger.warning('Detected another node type: '+str(node['type']))
                elif node['type']=='reaction':
                    self.logger.warning('Reation '+str(spe_name)+' contains species as predecessor')
                else:
                    self.logger.warning('Detected another node type: '+str(node['type']))
        return ordered_reac


    ##
    #
    #
    def orderedRetroReaction(self):
        end_species = self.endSpecies()
        if not len(end_species)==1:
            logging.warning('There are multiple end species: '+str(end_species))
        return self._countBackReac(end_species[0])



    ##
    #
    #
    def orderedReaction(self):
        start_species = self.startSpecies()
        if not len(start_species)==1:
            logging.warning('There are multiple end species: '+str(start_species))
        return self._countForwReac(start_species[0])



    ##
    #
    #
    def endSpecies(self):
        endspe = []
        self.logger.info(self._onlyProducedSpecies())
        for spe in self._onlyProducedSpecies():
            if len(self._countBackReac(spe))==self.num_reactions:
                endspe.append(spe)
        return endspe


    ##
    #
    #
    def startSpecies(self):
        staspe = []
        self.logger.info(self._onlyConsumedSpecies())
        for spe in self._onlyConsumedSpecies():
            if len(self._countForwReac(spe))==self.num_reactions:
                staspe.append(spe)
        return staspe
