import numpy as np
import tempfile
import logging
import pandas as pd
import rpSBML

logging.basicConfig(
    level=logging.DEBUG,
    #level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

class rpMerge:
    def __init__(self, 
                 species_group_id='central_species',
                 sink_species_group_id='rp_sink_species',
                 pathway_id='rp_pathway'):
        self.logger = logging.getLogger(__name__)
        self.pathway_id = pathway_id
        self.species_group_id = species_group_id
        self.sink_species_group_id = sink_species_group_id
        self.target_rpsbml = None
        self.source_rpsbml = None


    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################

    ## Check the libSBML calls
    #
    # Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html
    #
    # @param value The SBML call
    # @param message The string that describes the call
    def _checklibSBML(self, value, message):
        if value is None:
            self.logger.error('LibSBML returned a null value trying to ' + message + '.')
            raise AttributeError
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                self.logger.error(err_msg)
                raise AttributeError
        else:
            #self.logger.debug(message)
            return None


    ## Function to find the unique species
    #
    # pd_matrix is organised such that the rows are the simulated species and the columns are the measured ones
    #
    def _findUniqueRowColumn(self, pd_matrix):
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


    #######################################################################
    ###################################### INPUT FUNCTIONS ################
    #######################################################################




    ##########################################################################################
    #################################### REACTION ############################################
    ##########################################################################################


    ## Compare individual reactions and see if the measured pathway is contained within the simulated one
    #
    # Note that we assure that the match is 1:1 between species using the species match
    #
    def compareReaction(self, species_match, meas_reac, sim_reac):
        scores = []
        all_match = True
        ########### reactants #######
        ignore_reactants = []
        for meas_reactant in meas_reac.getListOfReactants():
            if meas_reactant.species in species_match:
                spe_found = False
                for sim_reactant in sim_reac.getListOfReactants():
                    if sim_reactant.species in species_match[meas_reactant.species] and not sim_reactant.species in ignore_reactants:
                        scores.append(species_match[meas_reactant.species][sim_reactant.species])
                        ignore_reactants.append(sim_reactant.species)
                        spe_found = True
                        break
                if not spe_found:
                    scores.append(0.0)
                    all_match = False
            else:
                logging.debug('Cannot find the measured species '+str(meas_reactant.species)+' in the the matched species: '+str(species_match))
                scores.append(0.0)
                all_match = False
        #products
        ignore_products = []
        for meas_product in meas_reac.getListOfProducts():
            if meas_product.species in species_match:
                pro_found = False
                for sim_product in sim_reac.getListOfProducts():
                    if sim_product.species in species_match[meas_product.species] and not sim_product.species in ignore_products:
                        scores.append(species_match[meas_product.species][sim_product.species])
                        ignore_products.append(sim_product.species)
                        pro_found = True
                        break
                if not pro_found:
                    scores.append(0.0)
                    all_match = False
            else:
                logging.debug('Cannot find the measured species '+str(meas_product.species)+' in the the matched species: '+str(species_match))
                scores.append(0.0)
                all_match = False
        return np.mean(scores), all_match


    ##########################################################################################
    ##################################### SPECIES ############################################
    ##########################################################################################

    ## Match all the measured chemical species to the simulated chemical species between two SBML 
    #
    # TODO: for all the measured species compare with the simualted one. Then find the measured and simulated species that match the best and exclude the 
    # simulated species from potentially matching with another
    #
    def compareSpecies(self, measured_rpsbml, sim_rpsbml, measured_comp_id=None, sim_comp_id=None):
        ############## compare species ###################
        meas_sim = {}
        sim_meas = {}
        species_match = {}
        for measured_species in measured_rpsbml.model.getListOfSpecies():
            #skip the species that are not in the right compartmennt, if specified
            if measured_comp_id and not measured_species.getCompartment()==measured_comp_id:
                continue
            self.logger.debug('--- Trying to match chemical species: '+str(measured_species.getId())+' ---')
            meas_sim[measured_species.getId()] = {}
            species_match[measured_species.getId()] = {}
            #species_match[measured_species.getId()] = {'id': None, 'score': 0.0, 'found': False}
            #TODO: need to exclude from the match if a simulated chemical species is already matched with a higher score to another measured species
            for sim_species in sim_rpsbml.model.getListOfSpecies():
                #skip the species that are not in the right compartmennt, if specified
                if sim_comp_id and not sim_species.getCompartment()==sim_comp_id:
                    continue
                meas_sim[measured_species.getId()][sim_species.getId()] = {'score': 0.0, 'found': False}
                if not sim_species.getId() in sim_meas:
                    sim_meas[sim_species.getId()] = {}
                sim_meas[sim_species.getId()][measured_species.getId()] = {'score': 0.0, 'found': False}
                measured_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(measured_species.getAnnotation())
                sim_rpsbml_brsynth_annot = sim_rpsbml.readBRSYNTHAnnotation(sim_species.getAnnotation())
                measured_miriam_annot = sim_rpsbml.readMIRIAMAnnotation(measured_species.getAnnotation())
                sim_miriam_annot = sim_rpsbml.readMIRIAMAnnotation(sim_species.getAnnotation())
                #### MIRIAM ####
                if sim_rpsbml.compareMIRIAMAnnotations(measured_species.getAnnotation(), sim_species.getAnnotation()):
                    self.logger.debug('--> Matched MIRIAM: '+str(sim_species.getId()))
                    meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.4
                    #meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2+0.2*jaccardMIRIAM(sim_miriam_annot, measured_miriam_annot)
                    meas_sim[measured_species.getId()][sim_species.getId()]['found'] = True
                ##### InChIKey ##########
                #find according to the inchikey -- allow partial matches
                #compare either inchikey in brsynth annnotation or MIRIAM annotation
                #NOTE: here we prioritise the BRSynth annotation inchikey over the MIRIAM
                measured_inchikey_split = None
                sim_rpsbml_inchikey_split = None
                if 'inchikey' in measured_brsynth_annot: 
                    measured_inchikey_split = measured_brsynth_annot['inchikey'].split('-')
                elif 'inchikey' in measured_miriam_annot:
                    if not len(measured_miriam_annot['inchikey'])==1:
                        #TODO: handle mutliple inchikey with mutliple compare and the highest comparison value kept
                        self.logger.warning('There are multiple inchikey values, taking the first one: '+str(measured_miriam_annot['inchikey']))
                    measured_inchikey_split = measured_miriam_annot['inchikey'][0].split('-')
                if 'inchikey' in sim_rpsbml_brsynth_annot:
                    sim_rpsbml_inchikey_split = sim_rpsbml_brsynth_annot['inchikey'].split('-')
                elif 'inchikey' in sim_miriam_annot:
                    if not len(sim_miriam_annot['inchikey'])==1:
                        #TODO: handle mutliple inchikey with mutliple compare and the highest comparison value kept
                        self.logger.warning('There are multiple inchikey values, taking the first one: '+str(sim_rpsbml_brsynth_annot['inchikey']))
                    sim_rpsbml_inchikey_split = sim_miriam_annot['inchikey'][0].split('-')
                if measured_inchikey_split and sim_rpsbml_inchikey_split:
                    if measured_inchikey_split[0]==sim_rpsbml_inchikey_split[0]:
                        self.logger.debug('Matched first layer InChIkey: ('+str(measured_inchikey_split)+' -- '+str(sim_rpsbml_inchikey_split)+')')
                        meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                        if measured_inchikey_split[1]==sim_rpsbml_inchikey_split[1]:
                            self.logger.debug('Matched second layer InChIkey: ('+str(measured_inchikey_split)+' -- '+str(sim_rpsbml_inchikey_split)+')')
                            meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                            meas_sim[measured_species.getId()][sim_species.getId()]['found'] = True
                            if measured_inchikey_split[2]==sim_rpsbml_inchikey_split[2]:
                                self.logger.debug('Matched third layer InChIkey: ('+str(measured_inchikey_split)+' -- '+str(sim_rpsbml_inchikey_split)+')')
                                meas_sim[measured_species.getId()][sim_species.getId()]['score'] += 0.2
                                meas_sim[measured_species.getId()][sim_species.getId()]['found'] = True
                sim_meas[sim_species.getId()][measured_species.getId()]['score'] = meas_sim[measured_species.getId()][sim_species.getId()]['score']
                sim_meas[sim_species.getId()][measured_species.getId()]['found'] = meas_sim[measured_species.getId()][sim_species.getId()]['found']
        #build the matrix to send
        meas_sim_mat = {}
        for i in meas_sim:
            meas_sim_mat[i] = {}
            for y in meas_sim[i]:
                meas_sim_mat[i][y] = meas_sim[i][y]['score']
        unique = self._findUniqueRowColumn(pd.DataFrame(meas_sim_mat))
        self.logger.debug('findUniqueRowColumn:')
        self.logger.debug(unique)
        for meas in meas_sim:
            if meas in unique:
                species_match[meas] = {}
                for unique_spe in unique[meas]:
                    species_match[meas][unique_spe] = round(meas_sim[meas][unique[meas][0]]['score'], 5)
            else:
                self.logger.warning('Cannot find a species match for the measured species: '+str(meas))
        self.logger.debug('#########################')
        self.logger.debug('species_match:')
        self.logger.debug(species_match)
        self.logger.debug('-----------------------')
        return species_match

    ######################################################################################################################
    ############################################### EC NUMBER ############################################################
    ######################################################################################################################

    def compareEC(meas_reac_miriam, sim_reac_miriam):
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
            logging.debug('Measured: ')
            logging.debug(measured_frac_ec)
            logging.debug('Simulated: ')
            logging.debug(sim_frac_ec)
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
            logging.warning('One of the two reactions does not have any EC entries.\nMeasured: '+str(meas_reac_miriam)+' \nSimulated: '+str(sim_reac_miriam))
            return 0.0



    #############################################################################################################
    ############################################ MERGE ##########################################################
    #############################################################################################################


    ## Merge two models species and reactions using the annotations to recognise the same species and reactions
    #
    # The source model has to have both the GROUPS and FBC packages enabled in its SBML. The course must have a groups
    #called rp_pathway. If not use the readSBML() function to create a model
    # We add the reactions and species from the rpsbml to the target_model
    # 
    # @param target_model input libsbml model object where we will add the reactions and species from self.model
    # @param pathway_id String default is rp_pathway, name of the pathway id of the groups object
    # @param addOrphanSpecies Boolean Default False
    # @param bilevel_obj Tuple of size 2 with the weights associated with the targetSink and GEM objective function
    #
    #TODO: add a confidence in the merge using the score in 
    def mergeModels(self,
                    source_rpsbml,
                    target_rpsbml):
        #target_rpsbml.model = target_document.getModel()
        #Find the ID's of the similar target_rpsbml.model species
        ################ MODEL FBC ########################
        if not target_rpsbml.model.isPackageEnabled('fbc'):
            self._checklibSBML(target_rpsbml.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        if not source_rpsbml.model.isPackageEnabled('fbc'):
            self._checklibSBML(source_rpsbml.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
        target_fbc = target_rpsbml.model.getPlugin('fbc')
        source_fbc = source_rpsbml.model.getPlugin('fbc')
        #note sure why one needs to set this as False
        self._checklibSBML(source_rpsbml.document.setPackageRequired('fbc', False), 'enabling FBC package')
        ################ UNITDEFINITIONS ######
        #return the list of unit definitions id's for the target to avoid overwritting
        #WARNING: this means that the original unit definitions will be prefered over the new one
        target_unitDefID = [i.getId() for i in target_rpsbml.model.getListOfUnitDefinitions()]
        for source_unitDef in source_rpsbml.model.getListOfUnitDefinitions():
            if not source_unitDef.getId() in target_unitDefID: #have to compare by ID since no annotation
                #create a new unitDef in the target
                target_unitDef = target_rpsbml.model.createUnitDefinition()
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
        compartment_conv = {}
        for source_compartment in source_rpsbml.model.getListOfCompartments():
            found = False
            target_ids = [i.getId() for i in target_rpsbml.model.getListOfCompartments()]
            source_annotation = source_compartment.getAnnotation()
            if not source_annotation:
                self.logger.warning('No annotation for the source of compartment '+str(source_compartment.getId()))
                continue
            #compare by MIRIAM first
            for target_compartment in target_rpsbml.model.getListOfCompartments():
                target_annotation = target_compartment.getAnnotation()
                if not target_annotation:
                    self.logger.warning('No annotation for the target of compartment: '+str(target_compartment.getId()))
                    continue
                if source_rpsbml.compareMIRIAMAnnotations(source_annotation, target_annotation):
                    found = True
                    compartment_conv[source_compartment.getId()] = target_compartment.getId() 
                    break
            if not found:
                #if the id is not found, see if the ids already exists
                if source_compartment.getId() in target_ids:
                    compartment_conv[source_compartment.getId()] = source_compartment.getId()
                    found = True
                #if there is not MIRIAM match and the id's differ then add it
                else:
                    target_compartment = target_rpsbml.model.createCompartment()
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
                    compartment_conv[target_compartment.getId()] = target_compartment.getId() 
        self.logging.debug('compartment_conv: '+str(compartment_conv))
        ################ PARAMETERS ###########
        #WARNING: here we compare by ID
        targetParametersID = [i.getId() for i in target_rpsbml.model.getListOfParameters()]
        for source_parameter in source_rpsbml.model.getListOfParameters():
            if not source_parameter.getId() in targetParametersID:
                target_parameter = target_rpsbml.model.createParameter()
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
        #Try to detect the current model species
        ### loop through all the source and the target equivalent compounds and return the matches
        '''
        output example: out['MNXM4__64__MNXC3'] = {'M_o2_c': 0.9575}
        note that it may occur that multiple matches occur, key is the first model entry and value is the second
        '''
        model_species_convert = self.compareSpecies(measured_rpsbml, sim_rpsbml, measured_comp_id=None, sim_comp_id=None)
        '''
        model_species_convert = {} #this is passed to be used for reaction ID conversion
        toAdd_model_species_ids = [i.getId() for i in source_rpsbml.model.getListOfSpecies()]
        model_species_ids = [i.getId() for i in source_rpsbml.model.getListOfSpecies()]
        targetModel_species_ids = [i.getId() for i in target_rpsbml.model.getListOfSpecies()]
        for model_species in source_rpsbml.model.getListOfSpecies():
            if model_species.getId() in targetModel_species_ids:
                toAdd_model_species_ids.remove(model_species.getId())
                continue
            model_species_annotation = model_species.getAnnotation()
            if not model_species_annotation:
                self.logger.error('Cannot find annotations for species: '+str(model_species.getId()))
                continue
            for targetModel_species in target_rpsbml.model.getListOfSpecies():
                targetModel_species_annotation = targetModel_species.getAnnotation()
                if not targetModel_species_annotation:
                    self.logger.warning('Cannot find annotations for species: '+str(targetModel_species.getId()))
                    continue
                ## here we compare perfect matches of MIRIAM annotations
                if source_rpsbml.compareMIRIAMAnnotations(model_species_annotation, targetModel_species_annotation):
                    #if we find a match from MIRIAM and not from ID we assume that it supersedes the ID and need to replace it in reactions
                    model_species_convert[model_species.getId()] = targetModel_species.getId()
                    toAdd_model_species_ids.remove(model_species.getId())
                    break
        '''
        self.logger.debug('model_species_convert: '+str(model_species_convert))
        #add the new species
        for model_species_id in toAdd_model_species_ids:
            model_species = source_rpsbml.model.getSpecies(model_species_id)
            if not model_species:
                self.logger.error('Cannot retreive model species: '+str(model_species_id))
            self._checklibSBML(model_species, 'fetching source species')
            targetModel_species = target_rpsbml.model.createSpecies()
            self._checklibSBML(targetModel_species, 'creating species')
            self._checklibSBML(targetModel_species.setMetaId(model_species.getMetaId()),
                    'setting target metaId')
            self._checklibSBML(targetModel_species.setId(model_species.getId()),
                    'setting target id')
            self._checklibSBML(targetModel_species.setCompartment(compartment_conv[model_species.getCompartment()]),
                    'setting target compartment')
            self._checklibSBML(targetModel_species.setInitialConcentration(
                model_species.getInitialConcentration()),
                    'setting target initial concentration')
            self._checklibSBML(targetModel_species.setBoundaryCondition(
                model_species.getBoundaryCondition()),
                    'setting target boundary concentration')
            self._checklibSBML(targetModel_species.setHasOnlySubstanceUnits(
                model_species.getHasOnlySubstanceUnits()),
                    'setting target has only substance units')
            self._checklibSBML(targetModel_species.setBoundaryCondition(
                model_species.getBoundaryCondition()),
                    'setting target boundary condition')
            self._checklibSBML(targetModel_species.setConstant(model_species.getConstant()),
                'setting target constant')
            self._checklibSBML(targetModel_species.setAnnotation(model_species.getAnnotation()),
                'setting target annotation')
        ################ REACTIONS ###################
        # Here we go by ID and check if it exists. If it does we need to check the species are the same
        # and if not, add it changing the ID #model_species_convert
        toAdd_model_reaction_ids = [i.getId() for i in source_rpsbml.model.getListOfReactions()]
        model_reaction_ids = [i.getId() for i in source_rpsbml.model.getListOfReactions()]
        targetModel_reaction_ids = [i.getId() for i in target_rpsbml.model.getListOfReactions()]
        reac_replace = {}
        for model_reaction in source_rpsbml.model.getListOfReactions():
            if model_reaction.getId() in targetModel_reaction_ids:
                toAdd_model_reaction_ids.remove(model_reaction.getId())
                continue
            model_reaction_speciesID = []
            for reactant in model_reaction.getListOfReactants():
                if reactant.species in model_species_convert:
                    model_reaction_speciesID.append(model_species_convert[reactant.species])
                else:
                    model_reaction_speciesID.append(reactant.species)
            model_reaction_productsID = []
            for product in model_reaction.getListOfProducts():
                if product.species in model_species_convert:
                    model_reaction_productsID.append(model_species_convert[product.species])
                else:
                    model_reaction_productsID.append(product.species)
            for targetModel_reaction in target_rpsbml.model.getListOfReactions():
                #loop through target model reactions
                targetModel_reaction_speciesID = [i.species for i in targetModel_reaction.getListOfReactants()]
                targetModel_reaction_productsID = [i.species for i in targetModel_reaction.getListOfProducts()]
                if not set(model_reaction_speciesID)-set(targetModel_reaction_speciesID) and not set(model_reaction_productsID)-set(targetModel_reaction_productsID):
                    self.logger.debug('The reactions species and products are the same: '+str(model_reaction.getId())+' --- '+str(targetModel_reaction.getId()))
                    #reac_replace[targetModel_reaction.getId()] = model_reaction.getId()
                    reac_replace[model_reaction.getId()] = targetModel_reaction.getId()
                    toAdd_model_reaction_ids.remove(model_reaction.getId())
                    continue
        #add the new reactions
        for model_reaction_id in toAdd_model_reaction_ids:
            model_reaction = source_rpsbml.model.getReaction(model_reaction_id)
            if not model_reaction:
                self.logger.error('Cannot retreive model reaction: '+str(model_reaction_id))
            self._checklibSBML(model_reaction, 'fetching source reaction')
            target_reaction = target_rpsbml.model.createReaction()
            self._checklibSBML(target_reaction, 'create reaction')
            target_fbc = target_reaction.getPlugin('fbc')
            self._checklibSBML(target_fbc, 'fetching target FBC package')
            source_fbc = model_reaction.getPlugin('fbc')
            self._checklibSBML(source_fbc, 'fetching source FBC package')
            source_upperFluxBound = source_fbc.getUpperFluxBound()
            self._checklibSBML(source_upperFluxBound, 'fetching upper flux bound')
            self._checklibSBML(target_fbc.setUpperFluxBound(source_upperFluxBound),
                    'setting upper flux bound')
            source_lowerFluxBound = source_fbc.getLowerFluxBound()
            self._checklibSBML(source_lowerFluxBound, 'fetching lower flux bound')
            self._checklibSBML(target_fbc.setLowerFluxBound(source_lowerFluxBound),
                    'setting lower flux bound')
            self._checklibSBML(target_reaction.setId(model_reaction.getId()), 'set reaction id')
            self._checklibSBML(target_reaction.setName(model_reaction.getName()), 'set name')
            self._checklibSBML(target_reaction.setSBOTerm(model_reaction.getSBOTerm()),
                    'setting the reaction system biology ontology (SBO)') #set as process
            #TODO: consider having the two parameters as input to the function
            self._checklibSBML(target_reaction.setReversible(model_reaction.getReversible()),
                    'set reaction reversibility flag')
            self._checklibSBML(target_reaction.setFast(model_reaction.getFast()),
                    'set reaction "fast" attribute')
            self._checklibSBML(target_reaction.setMetaId(model_reaction.getMetaId()), 'setting species meta_id')
            self._checklibSBML(target_reaction.setAnnotation(model_reaction.getAnnotation()),
                    'setting annotation for source reaction')
            #Reactants
            for model_reaction_speciesID in [i.species for i in model_reaction.getListOfReactants()]:
                target_reactant = target_reaction.createReactant()
                self._checklibSBML(target_reactant, 'create target reactant')
                if model_reaction_speciesID in model_species_convert:
                    self._checklibSBML(target_reactant.setSpecies(
                        model_species_convert[model_reaction_speciesID]), 'assign reactant species')
                else:
                    self._checklibSBML(target_reactant.setSpecies(model_reaction_speciesID),
                        'assign reactant species')
                source_reactant = model_reaction.getReactant(model_reaction_speciesID)
                self._checklibSBML(source_reactant, 'fetch source reactant')
                self._checklibSBML(target_reactant.setConstant(source_reactant.getConstant()),
                        'set "constant" on species '+str(source_reactant.getConstant()))
                self._checklibSBML(target_reactant.setStoichiometry(source_reactant.getStoichiometry()),
                        'set stoichiometry ('+str(source_reactant.getStoichiometry)+')')
            #Products
            for model_reaction_productID in [i.species for i in model_reaction.getListOfProducts()]:
                target_reactant = target_reaction.createProduct()
                self._checklibSBML(target_reactant, 'create target reactant')
                if model_reaction_productID in model_species_convert:
                    self._checklibSBML(target_reactant.setSpecies(
                        model_species_convert[model_reaction_productID]), 'assign reactant product')
                else:
                    self._checklibSBML(target_reactant.setSpecies(model_reaction_productID),
                        'assign reactant product')
                source_reactant = model_reaction.getProduct(model_reaction_productID)
                self._checklibSBML(source_reactant, 'fetch source reactant')
                self._checklibSBML(target_reactant.setConstant(source_reactant.getConstant()),
                        'set "constant" on product '+str(source_reactant.getConstant()))
                self._checklibSBML(target_reactant.setStoichiometry(source_reactant.getStoichiometry()),
                        'set stoichiometry ('+str(source_reactant.getStoichiometry)+')')
        #### GROUPS #####
        #TODO loop through the groups to add them
        if not target_rpsbml.model.isPackageEnabled('groups'):
            self._checklibSBML(target_rpsbml.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(source_rpsbml.document.setPackageRequired('groups', False), 'enabling groups package')
        source_groups = source_rpsbml.model.getPlugin('groups')
        self._checklibSBML(source_groups, 'fetching the source model groups')
        target_groups = target_rpsbml.model.getPlugin('groups')
        self._checklibSBML(target_groups, 'fetching the target model groups')
        self.logger.debug('model_species_convert: '+str(model_species_convert))
        self.logger.debug('reac_replace: '+str(reac_replace))
        #TODO: this will overwrite two groups of the same id, need to change
        for group in source_groups.getListOfGroups():
            #for all the species that need to be converted, replace the ones that are
            #if the group is the species group, replace the ones detected from model_species_convert
            if group.getId()==self.species_group_id or group.getId()==self.sink_species_group_id:
                for spe_conv in model_species_convert:
                    foundIt = False
                    for member in group.getListOfMembers():
                        if member.getIdRef()==spe_conv:
                            member.setIdRef(model_species_convert[spe_conv])
                            foundIt = True
                            break
            elif group.getId()==self.pathway_id:
                for reac_conv in reac_replace:
                    foundIt = False
                    for member in group.getListOfMembers():
                        if member.getIdRef()==reac_conv:
                            member.setIdRef(reac_replace[reac_conv])
                            foundIt = True
                            break
                    '''
                    if not foundIt:
                        self.logger.warning('Could not find '+str(spe_conv)+' to replace with '+str(model_species_convert[spe_conv]))
                    '''
            self._checklibSBML(target_groups.addGroup(group),
                    'copy the source groups to the target groups')
        ###### TITLES #####
        target_rpsbml.model.setId(target_rpsbml.model.getId()+'__'+source_rpsbml.model.getId())
        target_rpsbml.model.setName(target_rpsbml.model.getName()+' merged with '+source_rpsbml.model.getId())
        '''
        if fillOrphanSpecies==True:
            self.fillOrphan(target_rpsbml, self.pathway_id, compartment_id)
        '''
        return model_species_convert, reac_replace


