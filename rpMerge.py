from rpGraph import rpGraph
import logging

class rpMerge(rpGraph):
    """Class containing all the functions required to merge two SBML files together or two rpSBML objects
    """
    def __init__(self
                 model_name=None,
                 document=None,
                 path=None,
                 is_gem_sbml=False,
                 pathway_id='rp_pathway',
                 central_species_group_id='central_species',
                 sink_species_group_id='rp_sink_species_id'): 
        """Constructor of the class

        Automatically constructs the network when calling the construtor

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
        .. automethod:: _findUniqueRowColumn 
        """
        super().__init__(model_name, document, path, is_gem_sbml, pathway_id, central_species_group_id, sink_species_group_id))
        self.logger = logging.getLogger(__name__)


    ##################################################################################
    ############################# STATIC METHODS #####################################
    ##################################################################################


    @staticmethod
    def mergeSBMLFiles(path_source,
                       path_target,
                       path_merge,
                       is_source_gem=False,
                       is_target_gem=True,
                       del_sp_pro=False,
                       del_sp_react=True,
                       upper_flux_bound=999999.0,
                       lower_flux_bound=0.0,
                       compartment_id='MNXC3',
                       pathway_id='rp_pathway',
                       central_species_group_id='central_species',
                       sink_species_group_id='rp_sink_species_id'):
        """Static function that merges two SBML files together

        :param path_source: Path of the source SBML file
        :param path_target: Path of the target SBML file
        :param path_merge: Path of the output SBML file
        :param is_source_gem: Determine if the source is a GEM (default: False)
        :param is_target_gem: Determine if the targer is a GEM (default: True)
        :param del_sp_pro: Delete the single parent product (default: False)
        :param del_sp_react: Delete the single pareant reactant (default: True)
        :param upper_flux_bound: Upper flux FBC bound for the created reactions (default: 999999.0)
        :param lower_flux_bound: Lower flux FBC bound for the created reactions (default: 0.0)
        :param compartment_id: The compartment of the SBML file to add reactions (default: MNXC3)
        :param pathway_id: The heterologous pathway Groups id (default: rp_pathway)
        :param central_species_group_id: The central species Groups id (default: central_species)
        :param sink_species_group_id: The sink species Groups id (default: rp_sink_species_id)

        :type path_source: str
        :type path_target: str
        :type path_merge: str
        :type is_source_gem: bool
        :type is_target_gem: bool
        :type del_sp_pro: bool
        :type del_sp_react: bool
        :type upper_flux_bound: float
        :type lower_flux_bound: float
        :type compartment_id: str
        :type pathway_id: str
        :type central_species_group_id: str
        :type sink_species_group_id: str

        :return: Success or failure of the function
        :rtype: bool
        """
        if not os.path.exists(path_target):
            logging.error('Target SBML file is invalid: '+str(path_target))
            return False
        if not os.path.exists(path_source):
            logging.error('Target SBML file is invalid: '+str(path_source))
            return False
        source_rpsbml = rpMerge('source',
                                path=path_source,
                                is_gem_sbml=is_source_gem,
                                pathway_id=pathway_id,
                                central_species_group_id=central_species_group_id,
                                sink_species_group_id=sink_species_group_id)
        target_rpsbml = rpMerge('target',
                                path=path_target,
                                is_gem_sbml=is_target_gem,
                                pathway_id=pathway_id,
                                central_species_group_id=central_species_group_id,
                                sink_species_group_id=sink_species_group_id)
        source_rpsbml.mergeModels(target_rpsbml.model,
                                  del_sp_pro,
                                  del_sp_react,
                                  upper_flux_bound,
                                  lower_flux_bound,
                                  compartment_id,
                                  path_merge)
        return True


    ##################################################################################
    ############################# PRIVATE METHODS ####################################
    ##################################################################################


    def _findUniqueRowColumn(self, pd_matrix):
        """Private function that takes the matrix of similarity scores between the reactions or species of two models and finds the unqiue matches

        pd_matrix is organised such that the rows are the simulated species and the columns are the measured ones

        :param pd_matrix: Matrix of reactions or species of two models

        :type pd_matrix: np.array

        :return: Dictionary of matches
        :rtype: dict
        """
        self.logger.debug(pd_matrix)
        to_ret = {}
        ######################## filter by the global top values ################
        self.logger.debug('################ Filter best #############')
        #transform to np.array
        x = pd_matrix.values
        #resolve the rouding issues to find the max
        x = np.around(x, decimals=5)
        #first round involves finding the highest values and if found set to 0.0 the rows and columns (if unique)
        top = np.where(x==np.max(x))
        #as long as its unique keep looping
        if np.count_nonzero(x)==0:
            return to_ret
        while len(top[0])==1 and len(top[1])==1:
            if np.count_nonzero(x)==0:
                return to_ret
            pd_entry = pd_matrix.iloc[[top[0][0]],[top[1][0]]]
            row_name = str(pd_entry.index[0])
            col_name = str(pd_entry.columns[0])
            if col_name in to_ret:
                self.logger.debug('Overwriting (1): '+str(col_name))
                self.logger.debug(x)
            to_ret[col_name] = [row_name]
            #delete the rows and the columns 
            self.logger.debug('==================')
            self.logger.debug('Column: '+str(col_name))
            self.logger.debug('Row: '+str(row_name))
            pd_matrix.loc[:, col_name] = 0.0
            pd_matrix.loc[row_name, :] = 0.0
            x = pd_matrix.values
            x = np.around(x, decimals=5)
            top = np.where(x==np.max(x))
            self.logger.debug(pd_matrix)
            self.logger.debug(top)
            self.logger.debug('==================')
        #################### filter by columns (measured) top values ##############
        self.logger.debug('################ Filter by column best ############')
        x = pd_matrix.values
        x = np.around(x, decimals=5)
        if np.count_nonzero(x)==0:
            return to_ret
        reloop = True
        while reloop:
            if np.count_nonzero(x)==0:
                return to_ret
            reloop = False
            for col in range(len(x[0])):
                if np.count_nonzero(x[:,col])==0:
                    continue
                top_row = np.where(x[:,col]==np.max(x[:,col]))[0]
                if len(top_row)==1:
                    top_row = top_row[0]
                    #if top_row==0.0:
                    #    continue
                    #check to see if any other measured pathways have the same or larger score (accross)
                    row = list(x[top_row, :])
                    #remove current score consideration
                    row.pop(col)
                    if max(row)>=x[top_row, col]:
                        self.logger.warning('For col '+str(col)+' there are either better or equal values: '+str(row))
                        self.logger.warning(x)
                        continue
                    #if you perform any changes on the rows and columns, then you can perform the loop again
                    reloop = True
                    pd_entry = pd_matrix.iloc[[top_row],[col]]
                    self.logger.debug('==================')
                    row_name = pd_entry.index[0]
                    col_name = pd_entry.columns[0]
                    self.logger.debug('Column: '+str(col_name))
                    self.logger.debug('Row: '+str(row_name))
                    if col_name in to_ret:
                        self.logger.debug('Overwriting (2): '+str(col_name))
                        self.logger.debug(pd_matrix.values)
                    to_ret[col_name] = [row_name]
                    #delete the rows and the columns 
                    pd_matrix.loc[:, col_name] = 0.0
                    pd_matrix.loc[row_name, :] = 0.0
                    x = pd_matrix.values
                    x = np.around(x, decimals=5)
                    self.logger.debug(pd_matrix)
                    self.logger.debug('==================')
        ################## laslty if there are multiple values that are not 0.0 then account for that ######
        self.logger.debug('################# get the rest ##########')
        x = pd_matrix.values
        x = np.around(x, decimals=5)
        if np.count_nonzero(x)==0:
            return to_ret
        for col in range(len(x[0])):
            if not np.count_nonzero(x[:,col])==0:
                top_rows = np.where(x[:,col]==np.max(x[:,col]))[0]
                if len(top_rows)==1:
                    top_row = top_rows[0]
                    pd_entry = pd_matrix.iloc[[top_row],[col]]
                    row_name = pd_entry.index[0]
                    col_name = pd_entry.columns[0]
                    if not col_name in to_ret:
                        to_ret[col_name] = [row_name]
                    else:
                        self.logger.warning('At this point should never have only one: '+str(x[:,col]))
                        self.logger.warning(x)
                else:
                    for top_row in top_rows:
                        pd_entry = pd_matrix.iloc[[top_row],[col]]
                        row_name = pd_entry.index[0]
                        col_name = pd_entry.columns[0]
                        if not col_name in to_ret:
                            to_ret[col_name] = []
                        to_ret[col_name].append(row_name)
        self.logger.debug(pd_matrix)
        self.logger.debug('###################')
        return to_ret


    ##################################################################################
    ############################# PUBLIC METHODS #####################################
    ##################################################################################


    def checkSingleParent(self,
                          del_sp_pro=False,
                          del_sp_react=False,
                          upper_flux_bound=999999.0,
                          lower_flux_bound=0.0,
                          compartment_id='MNXM3')
        """Check if there are any single parent species in a heterologous pathways and if there are, either delete them or add reaction to complete the heterologous pathway

        :param del_sp_pro: Define if to delete the products or create reaction that consume it
        :param del_sp_react: Define if to delete the reactants or create reaction that produce it
        :param upper_flux_bound: The upper flux bounds unit definitions default when adding new reaction (Default: 999999.0)
        :param lower_flux_bound: The lower flux bounds unit definitions default when adding new reaction (Defaul: 0.0)
        :param compartment_id: The id of the model compartment

        :type del_sp_pro: bool
        :type del_sp_react: bool
        :type upper_flux_bound: float
        :type lower_flux_bound: float
        :type compartment_id: str

        :rtype: bool
        :return: Success of failure of the function
        """
        consumed_species_nid = self.onlyConsumedSpecies()
        produced_species_nid = self.onlyProducedSpecies()
        if del_sp_pro:
            for pro in produced_species_nid:
                self._checklibSBML(self.model.removeSpecies(pro), 'removing the following product species: '+str(pro))
        else:
            for pro in produced_species_nid:
                step = {'rule_id': None,
                        'left': {pro.split('__')[0]: 1},
                        'right': {},
                        'step': None,
                        'sub_step': None,
                        'path_id': None,
                        'transformation_id': None,
                        'rule_score': None,
                        'rule_ori_reac': None}
                #note that here the pathwats are passed as NOT being part of the heterologous pathways and 
                #thus will be ignored when/if we extract the rp_pathway from the full GEM model
                self.createReaction(pro+'__consumption',
                                    upper_flux_bound,
                                    lower_flux_bound,
                                    step,
                                    compartment_id)
        if del_sp_react:
            for react in consumed_species_nid:
                self._checklibSBML(self.model.removeSpecies(react), 'removing the following reactant species: '+str(react))
        else:
            for react in consumed_species_nid:
                step = {'rule_id': None,
                        'left': {},
                        'right': {react.split('__')[0]: 1},
                        'step': None,
                        'sub_step': None,
                        'path_id': None,
                        'transformation_id': None,
                        'rule_score': None,
                        'rule_ori_reac': None}
                #note that here the pathwats are passed as NOT being part of the heterologous pathways and 
                #thus will be ignored when/if we extract the rp_pathway from the full GEM model
                self.createReaction(react+'__production',
                                    upper_flux_bound,
                                    lower_flux_bound,
                                    step,
                                    compartment_id)
        return True


    ###################### REACTION ##############################


    # TODO: need to remove from the list reactions simulated reactions that have matched
    # TODO: Remove. This assumes that reactions can match multiple times, when in fact its impossible
    def compareReactions(self, species_match, target_rpsbml, source_rpsbml):
        """Compare the reactions of two SBML files

        Compare that all the measured species of a reactions are found within sim species to match with a reaction.
        We assume that there cannot be two reactions that have the same species and reactants. This is maintained by SBML

        :param species_match: The species match dictionary returned by compareSpecies()
        :param target_rpsbml: The target rpSBMl object
        :param source_rpsbml: The source rpSBML object

        :type species_match: dict
        :type target_rpsbml: rpSBML
        :type source_rpsbml: rpSBML

        :return: The dictionary of the reaction matches
        :rtype: dict
        """
        ############## compare the reactions #######################
        #construct sim reactions with species
        self.logger.debug('------ Comparing reactions --------')
        #match the reactants and products conversion to sim species
        tmp_reaction_match = {}
        source_target = {}
        target_source = {}
        for source_reaction in source_rpsbml.model.getListOfReactions():
            source_reaction_miriam = source_rpsbml.readMIRIAMAnnotation(source_reaction.getAnnotation())
            ################ construct the dict transforming the species #######
            source_target[source_reaction.getId()] = {}
            tmp_reaction_match[source_reaction.getId()] = {}
            for target_reaction in target_rpsbml.model.getListOfReactions():
                if not target_reaction.getId() in target_source:
                    target_source[target_reaction.getId()] = {}
                target_source[target_reaction.getId()][source_reaction.getId()] = {}
                source_target[source_reaction.getId()][target_reaction.getId()] = {}
                self.logger.debug('\t=========== '+str(target_reaction.getId())+' ==========')
                self.logger.debug('\t+++++++ Species match +++++++')
                tmp_reaction_match[source_reaction.getId()][target_reaction.getId()] = {'reactants': {},
                                                                             'reactants_score': 0.0,
                                                                             'products': {},
                                                                             'products_score': 0.0,
                                                                             'species_score': 0.0,
                                                                             'species_std': 0.0,
                                                                             'species_reaction': None,
                                                                             'ec_score': 0.0,
                                                                             'ec_reaction': None,
                                                                             'score': 0.0,
                                                                             'found': False}
                target_reaction = target_rpsbml.model.getReaction(target_reaction.getId())
                sim_reactants_id = [reactant.species for reactant in target_reaction.getListOfReactants()]
                sim_products_id = [product.species for product in target_reaction.getListOfProducts()]
                ############ species ############
                self.logger.debug('\tspecies_match: '+str(species_match))
                self.logger.debug('\tspecies_match: '+str(species_match.keys()))
                self.logger.debug('\tsim_reactants_id: '+str(sim_reactants_id))
                self.logger.debug('\tmeasured_reactants_id: '+str([i.species for i in source_reaction.getListOfReactants()]))
                self.logger.debug('\tsim_products_id: '+str(sim_products_id))
                self.logger.debug('\tmeasured_products_id: '+str([i.species for i in source_reaction.getListOfProducts()]))
                #ensure that the match is 1:1
                #1)Here we assume that a reaction cannot have twice the same species
                cannotBeSpecies = []
                #if there is a match then we loop again since removing it from the list of potential matches would be appropriate
                keep_going = True
                while keep_going:
                    self.logger.debug('\t\t----------------------------')
                    keep_going = False
                    for reactant in source_reaction.getListOfReactants():
                        self.logger.debug('\t\tReactant: '+str(reactant.species))
                        #if a species match has been found AND if such a match has been found
                        founReaIDs = [tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants'][i]['id'] for i in tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants'] if not tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants'][i]['id']==None]
                        self.logger.debug('\t\tfounReaIDs: '+str(founReaIDs))
                        if reactant.species and reactant.species in species_match and not list(species_match[reactant.species].keys())==[] and not reactant.species in founReaIDs:
                            best_spe = [k for k, v in sorted(species_match[reactant.species].items(), key=lambda item: item[1], reverse=True)][0]
                            tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants'][reactant.species] = {'id': best_spe, 'score': species_match[reactant.species][best_spe], 'found': True}
                            cannotBeSpecies.append(best_spe)
                        elif not reactant.species in tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants']:
                            self.logger.warning('\t\tCould not find the following measured reactant in the matched species: '+str(reactant.species))
                            tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants'][reactant.species] = {'id': None, 'score': 0.0, 'found': False}
                    for product in source_reaction.getListOfProducts():
                        self.logger.debug('\t\tProduct: '+str(product.species))
                        foundProIDs = [tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products'][i]['id'] for i in tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products'] if not tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products'][i]['id']==None]
                        self.logger.debug('\t\tfoundProIDs: '+str(foundProIDs))
                        if product.species and product.species in species_match and not list(species_match[product.species].keys())==[] and not product.species in foundProIDs:
                            best_spe = [k for k, v in sorted(species_match[product.species].items(), key=lambda item: item[1], reverse=True)][0]
                            tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants'][product.species] = {'id': best_spe, 'score': species_match[product.species][best_spe], 'found': True}
                            cannotBeSpecies.append(best_spe)
                        elif not product.species in tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products']:
                            self.logger.warning('\t\tCould not find the following measured product in the matched species: '+str(product.species))
                            tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products'][product.species] = {'id': None, 'score': 0.0, 'found': False}
                    self.logger.debug('\t\tcannotBeSpecies: '+str(cannotBeSpecies))
                reactants_score = [tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants'][i]['score'] for i in tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants']]
                reactants_found = [tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants'][i]['found'] for i in tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants']]
                tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['reactants_score'] = np.mean(reactants_score)
                products_score = [tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products'][i]['score'] for i in tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products']]
                products_found = [tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products'][i]['found'] for i in tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products']]
                tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['products_score'] = np.mean(products_score)
                ### calculate pathway species score
                tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['species_score'] = np.mean(reactants_score+products_score)
                tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['species_std'] = np.std(reactants_score+products_score)
                tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['species_reaction'] = target_reaction.getId()
                tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['found'] = all(reactants_found+products_found)
                tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['score'] = tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['species_score']
                target_source[target_reaction.getId()][source_reaction.getId()] = tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['score']
                source_target[source_reaction.getId()][target_reaction.getId()] = tmp_reaction_match[source_reaction.getId()][target_reaction.getId()]['score']
        ### matrix compare #####
        unique = self._findUniqueRowColumn(pd.DataFrame(source_target))
        self.logger.debug('findUniqueRowColumn')
        self.logger.debug(unique)
        reaction_match = {}
        for meas in source_target:
            reaction_match[meas] = {'id': None, 'score': 0.0, 'found': False}
            if meas in unique:
                if len(unique[meas])>1:
                    self.logger.debug('Multiple values may match, choosing the first arbitrarily: '+str(unique))
                reaction_match[meas]['id'] = unique[meas]
                reaction_match[meas]['score'] = round(tmp_reaction_match[meas][unique[meas][0]]['score'], 5)
                reaction_match[meas]['found'] = tmp_reaction_match[meas][unique[meas][0]]['found']
        #### compile a reaction score based on the ec and species scores
        self.logger.debug(tmp_reaction_match)
        self.logger.debug(reaction_match)
        self.logger.debug('-------------------------------')
        return reaction_match


    #TODO: change this with a flag so that all the reactants and products are the same
    def containedReaction(self, species_source_target, source_reaction, target_reaction):
        """Compare individual reactions and see if the source reaction is contained within the target one

        species_source_target: {'MNXM4__64__MNXC3': {'M_o2_c': 1.0}, 'MNXM10__64__MNXC3': {'M_nadh_c': 1.0}, 'CMPD_0000000003__64__MNXC3': {}, 'TARGET_0000000001__64__MNXC3': {}, 'MNXM188__64__MNXC3': {'M_anth_c': 1.0}, 'BC_32877__64__MNXC3': {'M_nh4_c': 0.8}, 'BC_32401__64__MNXC3': {'M_nad_c': 0.2}, 'BC_26705__64__MNXC3': {'M_h_c': 1.0}, 'BC_20662__64__MNXC3': {'M_co2_c': 1.0}}
        the first keys are the source compartment ids
        the second key is the source species id
        the value is the target species id
        Note that we assure that the match is 1:1 between species using the species match

        :param species_source_target: The comparison dictionary between the species of two SBML files
        :param source_reaction: The target reaction
        :param target_reaction: The source reaction

        :type species_source_target: dict
        :type source_reaction: libsbml.Reaction
        :type target_reaction: libsbml.Reaction

        :return: The score of the match and the dict of the match in that order
        :rtype: tuple
        """
        scores = []
        all_match = True
        ########### reactants #######
        ignore_reactants = []
        for source_reactant in source_reaction.getListOfReactants():
            if source_reactant.species in species_source_target:
                spe_found = False
                for target_reactiontant in target_reaction.getListOfReactants():
                    if target_reactiontant.species in species_source_target[source_reactant.species] and not target_reactiontant.species in ignore_reactants:
                        scores.append(species_source_target[source_reactant.species][target_reactiontant.species])
                        ignore_reactants.append(target_reactiontant.species)
                        spe_found = True
                        break
                if not spe_found:
                    scores.append(0.0)
                    all_match = False
            else:
                self.logger.debug('Cannot find the source species '+str(source_reactant.species)+' in the target species: '+str(species_source_target))
                scores.append(0.0)
                all_match = False
        #products
        ignore_products = []
        for source_product in source_reaction.getListOfProducts():
            if source_product.species in species_source_target:
                pro_found = False
                for sim_product in target_reaction.getListOfProducts():
                    if sim_product.species in species_source_target[source_product.species] and not sim_product.species in ignore_products:
                        scores.append(species_source_target[source_product.species][sim_product.species])
                        ignore_products.append(sim_product.species)
                        pro_found = True
                        break
                if not pro_found:
                    scores.append(0.0)
                    all_match = False
            else:
                self.logger.debug('Cannot find the measured species '+str(source_product.species)+' in the the matched species: '+str(species_source_target))
                scores.append(0.0)
                all_match = False
        return np.mean(scores), all_match


    #TODO: change this with a flag so that all the reactants and products are the same
    def compareReaction(self, species_source_target, source_reaction, target_reaction):
        """Compare two reactions and elect that they are the same if they have exactly the same reactants and products

        species_source_target: {'MNXM4__64__MNXC3': {'M_o2_c': 1.0}, 'MNXM10__64__MNXC3': {'M_nadh_c': 1.0}, 'CMPD_0000000003__64__MNXC3': {}, 'TARGET_0000000001__64__MNXC3': {}, 'MNXM188__64__MNXC3': {'M_anth_c': 1.0}, 'BC_32877__64__MNXC3': {'M_nh4_c': 0.8}, 'BC_32401__64__MNXC3': {'M_nad_c': 0.2}, 'BC_26705__64__MNXC3': {'M_h_c': 1.0}, 'BC_20662__64__MNXC3': {'M_co2_c': 1.0}}
        the first keys are the source compartment ids
        the second key is the source species id
        the value is the target species id
        Note that we assure that the match is 1:1 between species using the species match

        :param species_source_target: The comparison dictionary between the species of two SBML files
        :param source_reaction: The target reaction
        :param target_reaction: The source reaction

        :type species_source_target: dict
        :type source_reaction: libsbml.Reaction
        :type target_reaction: libsbml.Reaction

        :return: The score of the match and boolean if its a match or not
        :rtype: tuple
        """
        scores = []
        source_reactants = [i.species for i in source_reaction.getListOfReactants()]
        target_reactants = []
        for i in target_reaction.getListOfReactants():
            if i.species in species_source_target:
                if not species_source_target[i.species]=={}:
                    #WARNING: Taking the first one arbitrarely
                    conv_spe = [y for y in species_source_target[i.species]][0]
                    target_reactants.append(conv_spe)
                    scores.append(species_source_target[i.species][conv_spe])
                else:
                    target_reactants.append(i.species)
                    scores.append(1.0)
            else:
                target_reactants.append(i.species)
                scores.append(1.0)
        source_products = [i.species for i in source_reaction.getListOfProducts()]
        target_products = []
        for i in target_reaction.getListOfReactants():
            if i.species in species_source_target:
                if not species_source_target[i.species]=={}:
                    #WARNING: Taking the first one arbitrarely
                    conv_spe = [y for y in species_source_target[i.species]][0]
                    target_products.append(conv_spe)
                    scores.append(species_source_target[i.species][conv_spe])
                else:
                    target_products.append(i.species)
                    scores.append(1.0)
            else:
                target_products.append(i.species)
                scores.append(1.0)
        '''
        self.logger.debug('source_reactants: '+str(source_reactants))
        self.logger.debug('target_reactants: '+str(target_reactants))
        self.logger.debug('source_products: '+str(source_products))
        self.logger.debug('target_products: '+str(target_products))
        self.logger.debug(set(source_reactants)-set(target_reactants))
        self.logger.debug(set(source_products)-set(target_products))
        '''
        if not set(source_reactants)-set(target_reactants) and not set(source_products)-set(target_products):
            return np.mean(scores), True
        else:
            return np.mean(scores), False


    ########################### SPECIES ###############################


    # TODO: for all the measured species compare with the simualted one. Then find the measured and simulated species that match the best and exclude the 
    # simulated species from potentially matching with another
    def compareSpecies(self, comp_source_target, source_model, target_model):
        """Match all the measured chemical species to the simulated chemical species between two SBML

        :param comp_source_target: The comparison dictionary between the compartment of two SBML files
        :param source_rpsbml: The source rpSBML
        :param target_rpsbml: The target rpSBML

        :type species_source_target: dict
        :type source_rpsbml: rpSBML
        :type target_rpsbml: rpSBML

        :return: The compartment match dictionary
        :rtype: dict
        """
        ############## compare species ###################
        source_target = {}
        target_source = {}
        species_match = {}
        for source_species in source_model.getListOfSpecies():
            self.logger.debug('--- Trying to match chemical species: '+str(source_species.getId())+' ---')
            source_target[source_species.getId()] = {}
            species_match[source_species.getId()] = {}
            #species_match[source_species.getId()] = {'id': None, 'score': 0.0, 'found': False}
            #TODO: need to exclude from the match if a simulated chemical species is already matched with a higher score to another measured species
            for target_species in target_model.getListOfSpecies():
                #skip the species that are not in the same compartment as the source
                if not target_species.getCompartment()==comp_source_target[source_species.getCompartment()]:
                    continue
                source_target[source_species.getId()][target_species.getId()] = {'score': 0.0, 'found': False}
                if not target_species.getId() in target_source:
                    target_source[target_species.getId()] = {}
                target_source[target_species.getId()][source_species.getId()] = {'score': 0.0, 'found': False}
                source_brsynth_annot = self.readBRSYNTHAnnotation(source_species.getAnnotation())
                target_brsynth_annot = self.readBRSYNTHAnnotation(target_species.getAnnotation())
                source_miriam_annot = self.readMIRIAMAnnotation(source_species.getAnnotation())
                target_miriam_annot = self.readMIRIAMAnnotation(target_species.getAnnotation())
                #### MIRIAM ####
                if self.compareMIRIAMAnnotations(source_species.getAnnotation(), target_species.getAnnotation()):
                    self.logger.debug('--> Matched MIRIAM: '+str(target_species.getId()))
                    source_target[source_species.getId()][target_species.getId()]['score'] += 0.4
                    #source_target[source_species.getId()][target_species.getId()]['score'] += 0.2+0.2*jaccardMIRIAM(target_miriam_annot, source_miriam_annot)
                    source_target[source_species.getId()][target_species.getId()]['found'] = True
                ##### InChIKey ##########
                #find according to the inchikey -- allow partial matches
                #compare either inchikey in brsynth annnotation or MIRIAM annotation
                #NOTE: here we prioritise the BRSynth annotation inchikey over the MIRIAM
                source_inchikey_split = None
                target_inchikey_split = None
                if 'inchikey' in source_brsynth_annot:
                    source_inchikey_split = source_brsynth_annot['inchikey'].split('-')
                elif 'inchikey' in source_miriam_annot:
                    if not len(source_miriam_annot['inchikey'])==1:
                        #TODO: handle mutliple inchikey with mutliple compare and the highest comparison value kept
                        self.logger.warning('There are multiple inchikey values, taking the first one: '+str(source_miriam_annot['inchikey']))
                    source_inchikey_split = source_miriam_annot['inchikey'][0].split('-')
                if 'inchikey' in target_brsynth_annot:
                    target_inchikey_split = target_brsynth_annot['inchikey'].split('-')
                elif 'inchikey' in target_miriam_annot:
                    if not len(target_miriam_annot['inchikey'])==1:
                        #TODO: handle mutliple inchikey with mutliple compare and the highest comparison value kept
                        self.logger.warning('There are multiple inchikey values, taking the first one: '+str(target_brsynth_annot['inchikey']))
                    target_inchikey_split = target_miriam_annot['inchikey'][0].split('-')
                if source_inchikey_split and target_inchikey_split:
                    if source_inchikey_split[0]==target_inchikey_split[0]:
                        self.logger.debug('Matched first layer InChIkey: ('+str(source_inchikey_split)+' -- '+str(target_inchikey_split)+')')
                        source_target[source_species.getId()][target_species.getId()]['score'] += 0.2
                        if source_inchikey_split[1]==target_inchikey_split[1]:
                            self.logger.debug('Matched second layer InChIkey: ('+str(source_inchikey_split)+' -- '+str(target_inchikey_split)+')')
                            source_target[source_species.getId()][target_species.getId()]['score'] += 0.2
                            source_target[source_species.getId()][target_species.getId()]['found'] = True
                            if source_inchikey_split[2]==target_inchikey_split[2]:
                                self.logger.debug('Matched third layer InChIkey: ('+str(source_inchikey_split)+' -- '+str(target_inchikey_split)+')')
                                source_target[source_species.getId()][target_species.getId()]['score'] += 0.2
                                source_target[source_species.getId()][target_species.getId()]['found'] = True
                target_source[target_species.getId()][source_species.getId()]['score'] = source_target[source_species.getId()][target_species.getId()]['score']
                target_source[target_species.getId()][source_species.getId()]['found'] = source_target[source_species.getId()][target_species.getId()]['found']
        #build the matrix to send
        source_target_mat = {}
        for i in source_target:
            source_target_mat[i] = {}
            for y in source_target[i]:
                source_target_mat[i][y] = source_target[i][y]['score']
        unique = self._findUniqueRowColumn(pd.DataFrame(source_target_mat))
        self.logger.debug('findUniqueRowColumn:')
        self.logger.debug(unique)
        for meas in source_target:
            if meas in unique:
                species_match[meas] = {}
                for unique_spe in unique[meas]:
                    species_match[meas][unique_spe] = round(source_target[meas][unique[meas][0]]['score'], 5)
            else:
                self.logger.warning('Cannot find a species match for the measured species: '+str(meas))
        self.logger.debug('#########################')
        self.logger.debug('species_match:')
        self.logger.debug(species_match)
        self.logger.debug('-----------------------')
        return species_match


    ######################### EC NUMBER ####################################


    def compareEC(meas_reac_miriam, sim_reac_miriam):
        """Compare two MIRIAM annotations and find the similarity of their EC number

        :param meas_reac_miriam: The annotation object of the source
        :param sim_reac_miriam: The annotation object of the target

        :type meas_reac_miriam: libsbml.XMLNode
        :type sim_reac_miriam: libsbml.XMLNode

        :return: The match score
        :rtype: float
        """
        #Warning we only match a single reaction at a time -- assume that there cannot be more than one to match at a given time
        if 'ec-code' in meas_reac_miriam and 'ec-code' in sim_reac_miriam:
            measured_frac_ec = [[y for y in ec.split('.') if not y=='-'] for ec in meas_reac_miriam['ec-code']]
            sim_frac_ec = [[y for y in ec.split('.') if not y=='-'] for ec in sim_reac_miriam['ec-code']]
            #complete the ec numbers with None to be length of 4
            for i in range(len(measured_frac_ec)):
                for y in range(len(measured_frac_ec[i]), 4):
                    measured_frac_ec[i].append(None)
            for i in range(len(sim_frac_ec)):
                for y in range(len(sim_frac_ec[i]), 4):
                    sim_frac_ec[i].append(None)
            self.logger.debug('Measured: ')
            self.logger.debug(measured_frac_ec)
            self.logger.debug('Simulated: ')
            self.logger.debug(sim_frac_ec)
            best_ec_compare = {'meas_ec': [], 'sim_ec': [], 'score': 0.0, 'found': False}
            for ec_m in measured_frac_ec:
                for ec_s in sim_frac_ec:
                    tmp_score = 0.0
                    for i in range(4):
                        if not ec_m[i]==None and not ec_s[i]==None:
                            if ec_m[i]==ec_s[i]:
                                tmp_score += 0.25
                                if i==2:
                                    best_ec_compare['found'] = True
                            else:
                                break
                    if tmp_score>best_ec_compare['score']:
                        best_ec_compare['meas_ec'] = ec_m
                        best_ec_compare['sim_ec'] = ec_s
                        best_ec_compare['score'] = tmp_score
            return best_ec_compare['score']
        else:
            self.logger.warning('One of the two reactions does not have any EC entries.\nMeasured: '+str(meas_reac_miriam)+' \nSimulated: '+str(sim_reac_miriam))
            return 0.0

    
    ################################## MERGE #########################################


    #TODO: add a confidence in the merge using the score in 
    #TODO: seperate the different parts so that others may use it
    def mergeModels(self,
                    input_model,
                    del_sp_pro=False,
                    del_sp_react=True,
                    upper_flux_bound=999999.0,
                    lower_flux_bound=0.0,
                    compartment_id='MNXM3',
                    output_path=None):
        """Merge two models species and reactions using the annotations to recognise the same species and reactions

        The source model has to have both the GROUPS and FBC packages enabled in its SBML. The course must have a groups
        called rp_pathway. If not use the readSBML() function to create a model
        We add the reactions and species from the rpsbml to the target_model. When the output_path is specified, the function
        will output the model to that path and if left to the default None, the self model will be overwritten.

        :param input_model: The target model to add the source to
        :param del_sp_pro: Delete the single parent product (default: False)
        :param del_sp_react: Delete the single pareant reactant (default: True)
        :param upper_flux_bound: Upper flux FBC bound for the created reactions (default: 999999.0)
        :param lower_flux_bound: Lower flux FBC bound for the created reactions (default: 0.0)
        :param compartment_id: The compartment of the SBML file to add reactions (default: MNXC3)
        :param output_model: The path to the output file. If None the function overwrites

        :type input_model: Union[libsbml.Model, rpSBML, libsbml.Document, str]
        :type del_sp_pro: bool
        :type del_sp_react: bool
        :type upper_flux_bound: float
        :type lower_flux_bound: float
        :type compartment_id: str
        :type output_model: str

        :rtype: tuple
        :returns: tuple (species_source_target, reactions_source_target)
            - str species_source_target is dict
            - str reactions_source_target is dict
        Tuple of dict where the first entry`target model object is the species source to target conversion and the second is the reaction source to target conversion
        """
        if output_path and type(input_model)==libsbml.Model:
            self.logger.warning('Cannot pass a libsbml.Model object and excpect a physical file output... Setting output_path to None')
            output_path = None
        if type(input_model)==libsbml.Model:
            target_model = input_model
        elif type(input_model)==libsbml.SBMLDocument:
            target_model = input_model.model
        elif type(input_model)==rpSBML:
            target_model = input_model.model
        elif type(input_model)==str:
            if os.path.exists(input_model):
                rpsbml = rpSBML(model_name=self.model_name, path=input_model)
                target_model = rpsbml.model
            else:
                self.logger.error('The input model was detected to be a string and thus path, but the file does not seem to exists: '+str(input_model))
                return False
        else:
            self.logger.error('The input must be either a libsbml.SBMLDocument, libsbml.Model or the path to a model')
            return False
        #target_model = target_document.getModel()
        #Find the ID's of the similar target_model species
        ################ MODEL FBC ########################
        if not target_model.isPackageEnabled('fbc'):
            self._checklibSBML(target_model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        if not self.model.isPackageEnabled('fbc'):
            self._checklibSBML(self.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        target_fbc = target_model.getPlugin('fbc')
        source_fbc = self.model.getPlugin('fbc')
        #note sure why one needs to set this as False
        #self._checklibSBML(source_rpsbml.document.setPackageRequired('fbc', False), 'enabling FBC package')
        ################ UNITDEFINITIONS ######
        #return the list of unit definitions id's for the target to avoid overwritting
        #WARNING: this means that the original unit definitions will be prefered over the new one
        target_unitDefID = [i.getId() for i in target_model.getListOfUnitDefinitions()]
        for source_unitDef in self.model.getListOfUnitDefinitions():
            if not source_unitDef.getId() in target_unitDefID: #have to compare by ID since no annotation
                #create a new unitDef in the target
                target_unitDef = target_model.createUnitDefinition()
                self._checklibSBML(target_unitDef, 'fetching target unit definition')
                #copy unitDef info to the target
                self._checklibSBML(target_unitDef.setId(source_unitDef.getId()),
                    'setting target unit definition ID')
                self._checklibSBML(target_unitDef.setAnnotation(source_unitDef.getAnnotation()),
                    'setting target unit definition Annotation')
                for source_unit in source_unitDef.getListOfUnits():
                    #copy unit info to the target unitDef
                    target_unit = target_unitDef.createUnit()
                    self._checklibSBML(target_unit, 'creating target unit')
                    self._checklibSBML(target_unit.setKind(source_unit.getKind()),
                        'setting target unit kind')
                    self._checklibSBML(target_unit.setExponent(source_unit.getExponent()),
                        'setting target unit exponent')
                    self._checklibSBML(target_unit.setScale(source_unit.getScale()),
                        'setting target unit scale')
                    self._checklibSBML(target_unit.setMultiplier(source_unit.getMultiplier()),
                        'setting target unit multiplier')
                target_unitDefID.append(source_unitDef.getId()) #add to the list to make sure its not added twice
        ################ COMPARTMENTS ###############
        # Compare by MIRIAM annotations
        #Note that key is source and value is target conversion
        comp_source_target = {}
        for source_compartment in self.model.getListOfCompartments():
            found = False
            target_ids = [i.getId() for i in target_model.getListOfCompartments()]
            source_annotation = source_compartment.getAnnotation()
            if not source_annotation:
                self.logger.warning('No annotation for the source of compartment '+str(source_compartment.getId()))
                continue
            #compare by MIRIAM first
            for target_compartment in target_model.getListOfCompartments():
                target_annotation = target_compartment.getAnnotation()
                if not target_annotation:
                    self.logger.warning('No annotation for the target of compartment: '+str(target_compartment.getId()))
                    continue
                if self.compareMIRIAMAnnotations(source_annotation, target_annotation):
                    found = True
                    comp_source_target[source_compartment.getId()] = target_compartment.getId()
                    break
            if not found:
                #if the id is not found, see if the ids already exists
                if source_compartment.getId() in target_ids:
                    comp_source_target[source_compartment.getId()] = source_compartment.getId()
                    found = True
                #if there is not MIRIAM match and the id's differ then add it
                else:
                    target_compartment = target_model.createCompartment()
                    self._checklibSBML(target_compartment, 'Creating target compartment')
                    self._checklibSBML(target_compartment.setMetaId(source_compartment.getMetaId()),
                            'setting target metaId')
                    #make sure that the ID is different
                    if source_compartment.getId()==target_compartment.getId():
                        self._checklibSBML(target_compartment.setId(source_compartment.getId()+'_sourceModel'),
                                'setting target id')
                    else:
                        self._checklibSBML(target_compartment.setId(source_compartment.getId()),
                                'setting target id')
                    self._checklibSBML(target_compartment.setName(source_compartment.getName()),
                            'setting target name')
                    self._checklibSBML(target_compartment.setConstant(source_compartment.getConstant()),
                            'setting target constant')
                    self._checklibSBML(target_compartment.setAnnotation(source_compartment.getAnnotation()),
                            'setting target annotation')
                    self._checklibSBML(target_compartment.setSBOTerm(source_compartment.getSBOTerm()),
                            'setting target annotation')
                    comp_source_target[target_compartment.getId()] = target_compartment.getId()
        self.logger.debug('comp_source_target: '+str(comp_source_target))
        ################ PARAMETERS ###########
        #WARNING: here we compare by ID
        #TODO: need to improve and use the values that already exist in the SBML model
        targetParametersID = [i.getId() for i in target_model.getListOfParameters()]
        self.logger.debug('targetParametersID: '+str(targetParametersID)) #BUG: some of the id's are not detected and are almost overwritten (libSBML to the rescue)
        for source_parameter in self.model.getListOfParameters():
            if not source_parameter.getId() in targetParametersID:
                target_parameter = target_model.createParameter()
                self._checklibSBML(target_parameter, 'creating target parameter')
                self._checklibSBML(target_parameter.setId(source_parameter.getId()), 'setting target parameter ID')
                self._checklibSBML(target_parameter.setSBOTerm(source_parameter.getSBOTerm()),
                    'setting target parameter SBO')
                self._checklibSBML(target_parameter.setUnits(source_parameter.getUnits()),
                    'setting target parameter Units')
                self._checklibSBML(target_parameter.setValue(source_parameter.getValue()),
                    'setting target parameter Value')
                self._checklibSBML(target_parameter.setConstant(source_parameter.getConstant()),
                    'setting target parameter ID')
        ################ FBC GENE PRODUCTS ########################
        #WARNING: here we compare by ID
        targetGenProductID = [i.getId() for i in target_fbc.getListOfGeneProducts()]
        for source_geneProduct in source_fbc.getListOfGeneProducts():
            if not source_geneProduct.getId() in targetGenProductID:
                target_geneProduct = target_fbc.createGeneProduct()
                self._checklibSBML(target_geneProduct, 'creating target gene product')
                self._checklibSBML(target_geneProduct.setId(source_geneProduct.getId()),
                    'setting target gene product id')
                self._checklibSBML(target_geneProduct.setLabel(source_geneProduct.getLabel()),
                    'setting target gene product label')
                self._checklibSBML(target_geneProduct.setName(source_geneProduct.getName()),
                    'setting target gene product name')
                self._checklibSBML(target_geneProduct.setMetaId(source_geneProduct.getMetaId()),
                    'setting target gene product meta_id')
        ############### FBC OBJECTIVES ############
        #WARNING: here we compare by ID
        #TODO: if overlapping id's need to replace the id with modified, as for the species
        targetObjectiveID = [i.getId() for i in target_fbc.getListOfObjectives()]
        sourceObjectiveID = [i.getId() for i in source_fbc.getListOfObjectives()]
        for source_objective in source_fbc.getListOfObjectives():
            if not source_objective.getId() in targetObjectiveID:
                target_objective = target_fbc.createObjective()
                self._checklibSBML(target_objective, 'creating target objective')
                self._checklibSBML(target_objective.setId(source_objective.getId()), 'setting target objective')
                self._checklibSBML(target_objective.setName(source_objective.getName()), 'setting target objective')
                self._checklibSBML(target_objective.setType(source_objective.getType()),
                        'setting target objective type')
                for source_fluxObjective in source_objective.getListOfFluxObjectives():
                    target_fluxObjective = target_objective.createFluxObjective()
                    self._checklibSBML(target_fluxObjective, 'creating target flux objective')
                    self._checklibSBML(target_fluxObjective.setName(source_fluxObjective.getName()),
                        'setting target flux objective name')
                    self._checklibSBML(target_fluxObjective.setCoefficient(source_fluxObjective.getCoefficient()),
                        'setting target flux objective coefficient')
                    self._checklibSBML(target_fluxObjective.setReaction(source_fluxObjective.getReaction()),
                        'setting target flux objective reaction')
                    self._checklibSBML(target_fluxObjective.setAnnotation(source_fluxObjective.getAnnotation()),
                        'setting target flux obj annotation from source flux obj')
                self._checklibSBML(target_objective.setAnnotation(source_objective.getAnnotation()),
                        'setting target obj annotation from source obj')
        self.logger.debug('targetObjectiveID: '+str(targetObjectiveID))
        self.logger.debug('sourceObjectiveID: '+str(sourceObjectiveID))
        ################ SPECIES ####################
        species_source_target = self.compareSpecies(comp_source_target, self.model, target_model)
        self.logger.debug('species_source_target: '+str(species_source_target))
        target_species_ids = [i.id for i in target_model.getListOfSpecies()]
        for source_species in species_source_target:
            list_target = [i for i in species_source_target[source_species]]
            if source_species in list_target:
                self.logger.warning('The source ('+str(source_species)+') and target species ids ('+str(list_target)+') are the same')
            #if match, replace the annotation from the source to the target
            if not species_source_target[source_species]=={}:
                list_species = [i for i in species_source_target[source_species]]
                self.logger.debug('list_species: '+str(list_species))
                if len(list_species)==0:
                    continue
                    #self.logger.warning('Source species '+str(member.getIdRef())+' has been created in the target model')
                elif len(list_species)>1:
                    self.logger.warning('There are multiple matches to the species '+str(member.getIdRef())+'... taking the first one: '+str(list_species))
                #TODO: loop throught the annotations and replace the non-overlapping information
                target_member = target_model.getSpecies(list_species[0])
                source_member = self.model.getSpecies(source_species)
                self._checklibSBML(target_member, 'Retraiving the target species: '+str(list_species[0]))
                self._checklibSBML(source_member, 'Retreiving the source species: '+str(source_species))
                self._checklibSBML(target_member.setAnnotation(source_member.getAnnotation()), 'Replacing the annotations')
            #if no match then add it to the target model
            else:
                self.logger.debug('Creating source species '+str(source_species)+' in target rpsbml')
                source_species = self.model.getSpecies(source_species)
                if not source_species:
                    self.logger.error('Cannot retreive model species: '+str(source_species))
                else:
                    self._checklibSBML(source_species, 'fetching source species')
                    targetModel_species = target_model.createSpecies()
                    self._checklibSBML(targetModel_species, 'creating species')
                    self._checklibSBML(targetModel_species.setMetaId(source_species.getMetaId()),
                            'setting target metaId')
                    ## need to check if the id of the source species does not already exist in the target model
                    if source_species.getId() in target_species_ids:
                        target_species_id = self.model.id+'__'+str(source_species.getId())
                        if not source_species.getId() in species_source_target:
                            species_source_target[source_species.getId()] = {}
                        species_source_target[source_species.getId()][self.model.id+'__'+str(source_species.getId())] = 1.0
                    else:
                        target_species_id = source_species.getId()
                    self._checklibSBML(targetModel_species.setId(target_species_id),
                            'setting target id')
                    self._checklibSBML(targetModel_species.setCompartment(comp_source_target[source_species.getCompartment()]),
                            'setting target compartment')
                    self._checklibSBML(targetModel_species.setInitialConcentration(
                        source_species.getInitialConcentration()),
                            'setting target initial concentration')
                    self._checklibSBML(targetModel_species.setBoundaryCondition(
                        source_species.getBoundaryCondition()),
                            'setting target boundary concentration')
                    self._checklibSBML(targetModel_species.setHasOnlySubstanceUnits(
                        source_species.getHasOnlySubstanceUnits()),
                            'setting target has only substance units')
                    self._checklibSBML(targetModel_species.setBoundaryCondition(
                        source_species.getBoundaryCondition()),
                            'setting target boundary condition')
                    self._checklibSBML(targetModel_species.setConstant(source_species.getConstant()),
                        'setting target constant')
                    self._checklibSBML(targetModel_species.setAnnotation(source_species.getAnnotation()),
                        'setting target annotation')
        ################ REACTIONS ###################
        #TODO; consider the case where two reactions have the same ID's but are not the same reactions
        #TODO: if overlapping id's need to replace the id with modified, as for the species
        reactions_source_target = {}
        for source_reaction in self.model.getListOfReactions():
            is_found = False
            for target_reaction in target_model.getListOfReactions():
                score, match = self.compareReaction(species_source_target, source_reaction, target_reaction)
                if match:
                    self.logger.debug('Source reaction '+str(source_reaction)+' matches with target reaction '+str(target_reaction))
                    #source_reaction[source_reaction.getId()] = target_reaction.getId()
                    reactions_source_target[source_reaction.getId()] = target_reaction.getId()
                    is_found = True
                    break
            if not is_found:
                self.logger.debug('Cannot find source reaction: '+str(source_reaction.getId()))
                self._checklibSBML(source_reaction, 'fetching source reaction')
                target_reaction = target_model.createReaction()
                self._checklibSBML(target_reaction, 'create reaction')
                target_fbc = target_reaction.getPlugin('fbc')
                self._checklibSBML(target_fbc, 'fetching target FBC package')
                source_fbc = source_reaction.getPlugin('fbc')
                self._checklibSBML(source_fbc, 'fetching source FBC package')
                source_upperFluxBound = source_fbc.getUpperFluxBound()
                self._checklibSBML(source_upperFluxBound, 'fetching upper flux bound')
                self._checklibSBML(target_fbc.setUpperFluxBound(source_upperFluxBound),
                        'setting upper flux bound')
                source_lowerFluxBound = source_fbc.getLowerFluxBound()
                self._checklibSBML(source_lowerFluxBound, 'fetching lower flux bound')
                self._checklibSBML(target_fbc.setLowerFluxBound(source_lowerFluxBound),
                        'setting lower flux bound')
                self._checklibSBML(target_reaction.setId(source_reaction.getId()), 'set reaction id')
                self._checklibSBML(target_reaction.setName(source_reaction.getName()), 'set name')
                self._checklibSBML(target_reaction.setSBOTerm(source_reaction.getSBOTerm()),
                        'setting the reaction system biology ontology (SBO)') #set as process
                #TODO: consider having the two parameters as input to the function
                self._checklibSBML(target_reaction.setReversible(source_reaction.getReversible()),
                        'set reaction reversibility flag')
                self._checklibSBML(target_reaction.setFast(source_reaction.getFast()),
                        'set reaction "fast" attribute')
                self._checklibSBML(target_reaction.setMetaId(source_reaction.getMetaId()), 'setting species meta_id')
                self._checklibSBML(target_reaction.setAnnotation(source_reaction.getAnnotation()),
                        'setting annotation for source reaction')
                #Reactants
                self.logger.debug('Setting reactants')
                for source_reaction_reactantID in [i.species for i in source_reaction.getListOfReactants()]:
                    self.logger.debug('\tAdding '+str(source_reaction_reactantID))
                    target_reactant = target_reaction.createReactant()
                    self._checklibSBML(target_reactant, 'create target reactant')
                    if source_reaction_reactantID in species_source_target:
                        if not species_source_target[source_reaction_reactantID]=={}:
                            if len(species_source_target[source_reaction_reactantID])>1:
                                self.logger.warning('Multiple matches for '+str(source_reaction_reactantID)+': '+str(species_source_target[source_reaction_reactantID]))
                                self.logger.warning('Taking one the first one arbitrarely: '+str([i for i in species_source_target[source_reaction_reactantID]][0]))
                            #WARNING: taking the first one arbitrarely 
                            self._checklibSBML(target_reactant.setSpecies(
                                [i for i in species_source_target[source_reaction_reactantID]][0]), 'assign reactant species')
                        else:
                            self._checklibSBML(target_reactant.setSpecies(source_reaction_reactantID),
                                'assign reactant species')
                    else:
                        self._checklibSBML(target_reactant.setSpecies(source_reaction_reactantID),
                            'assign reactant species')
                    source_reactant = source_reaction.getReactant(source_reaction_reactantID)
                    self._checklibSBML(source_reactant, 'fetch source reactant')
                    self._checklibSBML(target_reactant.setConstant(source_reactant.getConstant()),
                            'set "constant" on species '+str(source_reactant.getConstant()))
                    self._checklibSBML(target_reactant.setStoichiometry(source_reactant.getStoichiometry()),
                            'set stoichiometry ('+str(source_reactant.getStoichiometry)+')')
                #Products
                self.logger.debug('Setting products')
                for source_reaction_productID in [i.species for i in source_reaction.getListOfProducts()]:
                    self.logger.debug('\tAdding '+str(source_reaction_productID))
                    target_product = target_reaction.createProduct()
                    self._checklibSBML(target_product, 'create target reactant')
                    if source_reaction_productID in species_source_target:
                        if not species_source_target[source_reaction_productID]=={}:
                            if len(species_source_target[source_reaction_reactantID])>1:
                                self.logger.warning('Multiple matches for '+str(source_reaction_productID)+': '+str(species_source_target[source_reaction_productID]))
                                self.logger.warning('Taking one arbitrarely')
                            #WARNING: taking the first one arbitrarely 
                            self._checklibSBML(target_product.setSpecies(
                                [i for i in species_source_target[source_reaction_productID]][0]), 'assign reactant product')
                        else:
                            self._checklibSBML(target_product.setSpecies(source_reaction_productID),
                                'assign reactant product')
                    else:
                        self._checklibSBML(target_product.setSpecies(source_reaction_productID),
                            'assign reactant product')
                    source_product = source_reaction.getProduct(source_reaction_productID)
                    self._checklibSBML(source_product, 'fetch source reactant')
                    self._checklibSBML(target_product.setConstant(source_product.getConstant()),
                            'set "constant" on product '+str(source_product.getConstant()))
                    self._checklibSBML(target_product.setStoichiometry(source_product.getStoichiometry()),
                            'set stoichiometry ('+str(source_product.getStoichiometry)+')')
        #### GROUPS #####
        #TODO loop through the groups to add them
        if not target_model.isPackageEnabled('groups'):
            self._checklibSBML(target_model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
        #!!!! must be set to false for no apparent reason
        #self._checklibSBML(source_rpsbml.document.setPackageRequired('groups', False), 'enabling groups package')
        source_groups = self.model.getPlugin('groups')
        self._checklibSBML(source_groups, 'fetching the source model groups')
        target_groups = target_model.getPlugin('groups')
        self._checklibSBML(target_groups, 'fetching the target model groups')
        #self.logger.debug('species_source_target: '+str(species_source_target))
        #self.logger.debug('reactions_source_target: '+str(reactions_source_target))
        source_groups_ids = [i.id for i in source_groups.getListOfGroups()]
        target_groups_ids = [i.id for i in target_groups.getListOfGroups()]
        #NOTE: only need to update the source species since these are the ones that are replaced with their equivalent
        for source_group in source_groups.getListOfGroups():
            #overwrite in the group the reaction members that have been replaced
            for member in source_group.getListOfMembers():
                if member.getIdRef() in reactions_source_target:
                    if reactions_source_target[member.getIdRef()]:
                        member.setIdRef(reactions_source_target[member.getIdRef()])
            #overwrite in the group the species members that have been replaced
            for member in source_group.getListOfMembers():
                if member.getIdRef() in species_source_target:
                    if species_source_target[member.getIdRef()]:
                        list_species = [i for i in species_source_target[member.getIdRef()]]
                        self.logger.debug('species_source_target: '+str(species_source_target))
                        self.logger.debug('list_species: '+str(list_species))
                        if len(list_species)==0:
                            continue
                            #self.logger.warning('Source species '+str(member.getIdRef())+' has been created in the target model')
                        elif len(list_species)>1:
                            self.logger.warning('There are multiple matches to the species '+str(member.getIdRef())+'... taking the first one: '+str(list_species))
                        self._checklibSBML(member.setIdRef(list_species[0]), 'Setting name to the groups member')
            #create and add the groups if a source group does not exist in the target
            if not source_group.id in target_groups_ids:
                self._checklibSBML(target_groups.addGroup(source_group),
                    'copy the source groups to the target groups')
            #if the group already exists in the target then need to add new members
            else:
                target_group = target_groups.getGroup(source_group.id)
                target_group_ids = [i.getIdRef() for i in target_group.getListOfMembers()]
                for member in source_group.getListOfMembers():
                    if member.getIdRef() not in target_group_ids:
                        new_member = target_group.createMember()
                        self._checklibSBML(new_member, 'Creating a new groups member')
                        self._checklibSBML(new_member.setIdRef(member.getIdRef()), 'Setting name to the groups member')
        """
        for group in source_groups.getListOfGroups():
            #for all the species that need to be converted, replace the ones that are
            #if the group is the species group, replace the ones detected from species_source_target
            if group.getId()==species_group_id or group.getId()==sink_species_group_id:
                for member in group.getListOfMembers():
                    if member.getIdRef() in species_source_target:
                        list_species = [i for i in species_source_target[member.getIdRef()]]
                        self.logger.debug('species_source_target: '+str(species_source_target))
                        self.logger.debug('list_species: '+str(list_species))
                        if len(list_species)==0:
                            self.logger.warning('Source species '+str(member.getIdRef())+' has been created in the target model')
                        elif len(list_species)>1:
                            self.logger.warning('There are multiple matches to the species '+str(member.getIdRef())+'... taking the first one: '+str(list_species))
                            #WARNING: taking the first one arbitrarely 
                            member.setIdRef(list_species[0])
                        else:
                            member.setIdRef(list_species[0])
            elif group.getId()==pathway_id:
                for member in group.getListOfMembers():
                    if member.getIdRef() in reactions_source_target:
                        member.setIdRef(reactions_source_target[member.getIdRef()])
            self._checklibSBML(target_groups.addGroup(group),
                    'copy the source groups to the target groups')
        """
        ###### TITLES #####
        target_model.setId(target_model.getId()+'__'+self.model.getId())
        target_model.setName(target_model.getName()+' merged with '+self.model.getId())
        #detect single parent species and deal with them
		self.checkSingleParent(del_sp_pro, del_sp_react, upper_flux_bound, lower_flux_bound, compartment_id)
        ########### OUTPUT ##############
        #if the output is specified, then generate the different 
        if output_path:
            if type(input_model)==libsbml.SBMLDocument:
                libsbml.writeSBMLToFile(input_model, output_path)
            elif type(input_model)==rpSBML:
                input_model.writeSBML(output_path)
            elif type(input_model)==str:
                libsbml.writeSBMLToFile(rpsbml.document, output_path)
        else: #if there are no output specified, we overwrite the current model
            if type(input_model)==libsbml.SBMLDocument:
                self.document = rpsbml.document
                self.model = target_model
            elif type(input_model)==libsbml.Model:
                self.document = None
            elif type(input_model)==rpSBML:
                self.document = input_model.document
                self.model = input_model.model
                self.model_name = input_model.model_name
            elif type(input_model)==str:
                self.document = rpsbml.document
                self.model = rpsbml.model
        return species_source_target, reactions_source_target
