import libsbml
import tempfile
import glob
import logging


from .rpMerge import rpMerge


__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = []
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Melchior du Lac"
__status__ = "Development"

#logging.root.setLevel(logging.NOTSET)

logging.basicConfig(
    #level=logging.DEBUG,
    level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


#TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes
class rpFBA(rpMerge):
    import cobra
    from cobra.flux_analysis import pfba
    """Class to simulate an rpsbml object using different FBA types and objective functions
    """
    def __init__(self,
                 gem_sbml_path=None,
                 rpsbml_path=None,
                 del_sp_pro=False,
                 del_sp_react=True,
                 upper_flux_bound=999999.0,
                 lower_flux_bound=0.0,
                 compartment_id='MNXC3',
                 is_gem_sbml=False,
                 model_name=None,
                 pathway_id='rp_pathway',
                 rpcache=None):
        """Default constructor

        The class inherits from rpMerge. The unique parameter of this class is gem_sbml_path, that determines the path to the GEM model file
        that the internal SBML will merge to (controlled by path). The uppser and lower fluxes and compartment_id control the creation of novel reactions
        when faced with single parent species that we will fill. Note that if both del_sp_pro and del_sp_react are False, then these are ignored since
        single parent species are simply deleted. The different group_id are specified due to the rpGraph package (not that important). The rest of the
        input parameters relate to rpSBML (model_name, document, path....).

        NOTE: that is_gem_sbml is specified for the CURRENT model and not the passed one through gem_sbml_path

        .. document private functions
        .. automethod:: _convertToCobra
        .. automethod:: _writeAnalysisResults
        """
        super().__init__(is_gem_sbml=True,
                         model_name=model_name,
                         path=gem_sbml_path,
                         rpcache=rpcache)
        self.logger = logging.getLogger(__name__)
        self.logger.debug('Started instance of rpFBA')
        self.species_source_target = None
        self.reactions_source_target = None
        self.rpsbml_path = rpsbml_path
        self.gem_sbml_path = gem_sbml_path
        self.del_sp_pro = del_sp_pro
        self.del_sp_react = del_sp_react
        self.upper_flux_bound = upper_flux_bound
        self.lower_flux_bound = lower_flux_bound
        self.compartment_id = compartment_id
        self.is_gem_sbml = is_gem_sbml
        self.model_name = model_name
        self.pathway_id = pathway_id
        self.rpcache = rpcache
        self.cobra_model = None
        #TODO enable FBC if not done so
        if self.model:
            if gem_sbml_path:
                status = self.mergeModels(self.rpsbml_path,
                                          self.del_sp_pro,
                                          self.del_sp_react,
                                          self.upper_flux_bound,
                                          self.lower_flux_bound,
                                          self.compartment_id,
                                          self.pathway_id)
                if not status:
                    self.logger.warning('Failed to merge the models')


    ##########################################################
    ################# Private Functions ######################
    ##########################################################


    def _convertToCobra(self):
        """Convert the rpSBML object to cobra object

        :return: Success or failure of the function
        :rtype: bool
        """
        cobra_model = None
        try:
            with tempfile.TemporaryDirectory() as tmp_output_folder:
                self.writeSBML(tmp_output_folder)
                #logging.info(glob.glob(tmp_output_folder+'/*'))
                #logging.info(cobra.io.validate_sbml_model(glob.glob(tmp_output_folder+'/*')[0]))
                cobra_model = cobra.io.read_sbml_model(glob.glob(tmp_output_folder+'/*')[0], use_fbc_package=True)
            #self.cobra_model = cobra.io.read_sbml_model(self.document.toXMLNode().toXMLString(), use_fbc_package=True)
            #use CPLEX
            # self.cobra_model.solver = 'cplex'
        except cobra.io.sbml.CobraSBMLError as e:
            self.logger.error(e)
            self.logger.error('Cannot convert the libSBML model to Cobra')
            return False
        return cobra_model


    def _writeAnalysisResults(self, objective_id, cobra_results, pathway_id='rp_pathway'):
        """Method to harcode into BRSynth annotations the results of a COBRA analysis

        :param objective_id: The id of the objective to optimise
        :param cobra_results: The cobrapy results object
        :param pathway_id: The id of the heterologous pathway group (Default: rp_pathway)

        :type cobra_results: cobra.ModelSummary
        :type objective_id: str
        :type pathway_id: str

        :return: None
        :rtype: None
        """
        self.logger.debug('----- Setting the results for '+str(objective_id)+ ' -----')
        groups = self.model.getPlugin('groups')
        self._checklibSBML(groups, 'Getting groups plugin')
        try:
            rp_pathway = groups.getGroup(pathway_id)
            self._checklibSBML(rp_pathway, 'retreiving the group: '+str(pathway_id))
        except AttributeError:
            self.logger.warning('The group '+str(pathway_id)+' does not exist... creating it')
            self.createPathway(pathway_id)
            rp_pathway = groups.getGroup(pathway_id)
            self._checklibSBML(rp_pathway, 'Getting RP pathway')
        #write the results to the rp_pathway
        self.logger.debug('Set '+str(pathway_id)+' with '+str('fba_'+str(objective_id))+' to '+str(cobra_results.objective_value))
        self.addUpdateBRSynth(rp_pathway, 'fba_'+str(objective_id), str(cobra_results.objective_value), 'mmol_per_gDW_per_hr', False)
        #get the objective
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC plugin')
        obj = fbc_plugin.getObjective(objective_id)
        self._checklibSBML(obj, 'Getting objective '+str(objective_id))
        self.addUpdateBRSynth(obj, 'flux_value', str(cobra_results.objective_value), 'mmol_per_gDW_per_hr', False)
        self.logger.debug('Set the objective '+str(objective_id)+' a flux_value of '+str(cobra_results.objective_value))
        for flux_obj in obj.getListOfFluxObjectives():
            #sometimes flux cannot be returned
            if cobra_results.fluxes.get(flux_obj.getReaction())==None:
                self.logger.warning('Cobra BUG: Cannot retreive '+str(flux_obj.getReaction())+' flux from cobrapy... setting to 0.0')
                self.addUpdateBRSynth(flux_obj, 'flux_value', str(0.0), 'mmol_per_gDW_per_hr', False)
                self.logger.debug('Set the reaction '+str(flux_obj.getReaction())+' a flux_value of '+str(0.0))
            else:
                self.addUpdateBRSynth(flux_obj, 'flux_value', str(cobra_results.fluxes.get(flux_obj.getReaction())), 'mmol_per_gDW_per_hr', False)
                self.logger.debug('Set the reaction '+str(flux_obj.getReaction())+' a flux_value of '+str(cobra_results.fluxes.get(flux_obj.getReaction())))
        #write all the results to the reactions of pathway_id
        for member in rp_pathway.getListOfMembers():
            try:
                self._checklibSBML(self.model.getReaction(member.getIdRef()), 'retreiving the reaction: '+str(member.getIdRef()))
            except AttributeError:
                self.logger.error('Cannot retreive the following reaction: '+str(member.getIdRef()))
                continue
            self.logger.debug('Set the reaction '+str(member.getIdRef())+' a '+str('fba_'+str(objective_id))+' of '+str(cobra_results.fluxes.get(member.getId())))
            self.addUpdateBRSynth(member, 'fba_'+str(objective_id), str(cobra_results.fluxes.get(member.getId())), 'mmol_per_gDW_per_hr', False)


    ##################################################################
    ######################## PUBLIC FUNCTIONS ########################
    ##################################################################


    def mergeModels(self,
                    rpsbml_path=None,
                    del_sp_pro=False,
                    del_sp_react=True,
                    upper_flux_bound=999999.0,
                    lower_flux_bound=0.0,
                    compartment_id='MNXC3',
                    pathway_id='rp_pathway'):
        if not del_sp_pro==False and not del_sp_pro:
            if not self.del_sp_pro==False and not self.del_sp_pro:
                self.logger.error('Must define either self ('+str(self.del_sp_pro)+') or input del_sp_pro ('+str(del_sp_pro)+')')
                return False
            d_s_p = self.del_sp_pro
        else:
            d_s_p = del_sp_pro
        if not del_sp_react==False and not del_sp_react:
            if not self.del_sp_react==False and not self.del_sp_react:
                self.logger.error('Must define either self ('+str(self.del_sp_react)+') or input del_sp_react ('+str(del_sp_react)+')')
                return False
            d_s_r = self.del_sp_react
        else:
            d_s_r = del_sp_react
        if not upper_flux_bound==0.0 and not upper_flux_bound:
            if not self.upper_flux_bound==0.0 and not self.upper_flux_bound:
                self.logger.error('Must define either self ('+str(self.upper_flux_bound)+') or input upper_flux_bound ('+str(upper_flux_bound)+')')
                return False
            u_f_b = self.upper_flux_bound
        else:
            u_f_b = upper_flux_bound
        if not lower_flux_bound==0.0 and not lower_flux_bound:
            if not self.lower_flux_bound==0.0 and not self.lower_flux_bound:
                self.logger.error('Must define either self ('+str(self.lower_flux_bound)+') or input lower_flux_bound ('+str(lower_flux_bound)+')')
                return False
            l_f_b = self.lower_flux_bound
        else:
            l_f_b = lower_flux_bound
        if not compartment_id:
            if not self.compartment_id:
                self.logger.error('Must define either self ('+str(self.compartment_id)+') or input compartment_id ('+str(compartment_id)+')')
                return False
            c_i = self.compartment_id
        else:
            c_i = compartment_id
        if not pathway_id:
            if not self.pathway_id:
                self.logger.error('Must define either self ('+str(self.pathway_id)+') or input pathway_id ('+str(pathway_id)+')')
                return False
            p_i = self.pathway_id
        else:
            p_i = pathway_id
        self.species_source_target, self.reactions_source_target = super().mergeModels(rpsbml_path, d_s_p, d_s_r, u_f_b, l_f_b, c_i, p_i)
        return True


    def writeSBML(self, path=None, keep_merged=True):
        """Export the metabolic network to a SBML filem either the fully merged versin or the rpSBML only version

        :param path: Path to the output SBML file

        :type path: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: bool
        :return: Success or failure of the command
        """
        ####### check the path #########
        #need to determine where are the path id's coming from
        if not keep_merged:
            if not self.path:
                self.logger.error('Must have path initiated to return the non-merged version of the model')
                return False
            #open the original file once again
            skinny_rpmerge = rpMerge(is_gem_sbml=False,
                                     path=self.rpsbml_path,
                                     rpcache=self.rpcache)
            #TODO: need to add the newly created reaction from the merge to the skinny rpSBML version
            #overwrite the original BRSynth annotations with the new ones
            skinny_rpmerge.updateBRSynthPathway(self.asDict())
            #WARNING: This will not keep the created reactions at the gafilling single parent process
            #add the objectives and overwtie as well
            source_fbc = self.model.getPlugin('fbc')
            self._checklibSBML(source_fbc, 'Getting source FBC')
            target_fbc = skinny_rpmerge.model.getPlugin('fbc')
            self._checklibSBML(target_fbc, 'Getting target FBC')
            target_objID = [i.getId() for i in target_fbc.getListOfObjectives()]
            for source_obj in source_fbc.getListOfObjectives():
                source_obj_id = source_obj.getId()
                if source_obj.getId() in target_objID:
                    target_obj = target_fbc.getObjective(source_obj.getId())
                    self._checklibSBML(target_obj, 'Getting target objective')
                    self._checklibSBML(target_obj.setAnnotation(source_obj.getAnnotation()), 'setting annotation')
                    for target_fluxObj in target_obj.getListOfFluxObjectives():
                        for source_fluxObj in source_obj.getListOfFluxObjectives():
                            if target_fluxObj.getReaction()==source_fluxObj.getReaction():
                                self._checklibSBML(target_fluxObj.setAnnotation(source_fluxObj.getAnnotation()), 'setting annotation')
                else:
                    self._checklibSBML(target_fbc.addObjective(source_obj), 'Adding objective')
            #rpsbml.createMultiFluxObj('obj_RP1_sink', ['RP1_sink'], [1])
            target_fbc.setActiveObjectiveId(source_obj_id) #tmp random assigenement of objective
            return skinny_rpmerge.writeSBML(path)
        else:
            return super().writeSBML(path)


    def runMultiObjective(self,
                          reactions,
                          coefficients,
                          is_max=True,
                          pathway_id='rp_pathway',
                          objective_id=None):
        """Run FBA using multiple objectives

        :param reactions: The ids of the reactions involved in the objective
        :param coefficients: The coefficients associated with the reactions id
        :param is_max: Maximise or minimise the objective (Default: True)
        :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
        :param objective_id: Overwrite the default id (Default: None)

        :type reactions: list
        :type coefficients: list
        :type is_max: bool
        :type pathway_id: str
        :type objective_id: str

        :return: Success or failure of the function
        :rtype: bool
        """
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        objective_id = self.findCreateObjective(reactions, coefficients, is_max)
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        if not self.cobra_model:
            self.cobra_model = self._convertToCobra()
        if not self.cobra_model:
            self._writeAnalysisResults(objective_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = self.cobra_model.optimize()
        self._writeAnalysisResults(objective_id, cobra_results, pathway_id)
        return 0.0, True


    #TODO
    #def runMultiObjectiveParsimonou


    def runFBA(self, reaction_id, coefficient=1.0, is_max=True, pathway_id='rp_pathway', objective_id=None):
        """Run FBA using a single objective

        :param reaction_id: The id of the reactions involved in the objective
        :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
        :param is_max: Maximise or minimise the objective (Default: True)
        :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
        :param objective_id: Overwrite the default id (Default: None)

        :type reaction_id: str
        :type coefficient: float
        :type is_max: bool
        :type pathway_id: str
        :type objective_id: str

        :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
        :rtype: tuple
        """
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        objective_id = self.findCreateObjective([reaction_id], [coefficient], is_max, objective_id)
        #run the FBA
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        if not self.cobra_model:
            self.cobra_model = self._convertToCobra()
        if not self.cobra_model:
            self._writeAnalysisResults(objective_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = self.cobra_model.optimize()
        self._writeAnalysisResults(objective_id, cobra_results, pathway_id)
        return cobra_results.objective_value, True


    def runParsimoniousFBA(self, reaction_id, coefficient=1.0, fraction_of_optimum=0.95, is_max=True, pathway_id='rp_pathway', objective_id=None):
        """Run parsimonious FBA using a single objective

        :param reaction_id: The id of the reactions involved in the objective
        :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
        :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.95)
        :param is_max: Maximise or minimise the objective (Default: True)
        :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
        :param objective_id: Overwrite the default id (Default: None)

        :type reaction_id: str
        :type coefficient: float
        :type fraction_of_optimum: float
        :type is_max: bool
        :type pathway_id: str
        :type objective_id: str

        :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
        :rtype: tuple
        """
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        objective_id = self.findCreateObjective([reaction_id], [coefficient], is_max, objective_id)
        #run the FBA
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        if not self.cobra_model:
            self.cobra_model = self._convertToCobra()
        if not self.cobra_model:
            self._writeAnalysisResults(objective_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = pfba(self.cobra_model, fraction_of_optimum)
        self._writeAnalysisResults(objective_id, cobra_results, pathway_id)
        return cobra_results.objective_value, True


    def runFractionReaction(self,
                            source_reaction,
                            source_coefficient,
                            target_reaction,
                            target_coefficient,
                            fraction_of_source=0.75,
                            is_max=True,
                            pathway_id='rp_pathway',
                            objective_id='obj_fraction'):
        """Optimise for a target reaction while fixing a source reaction to the fraction of its optimum

        :param source_reaction: The id of the source reaction
        :param source_coefficient: The source coefficient associated with the source reaction id
        :param target_reaction: The id of the target reaction
        :param target_coefficient: The source coefficient associated with the target reaction id
        :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.75)
        :param is_max: Maximise or minimise the objective (Default: True)
        :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
        :param objective_id: Overwrite the default id (Default: None)

        :type source_reaction: str
        :type source_coefficient: float
        :type target_reaction: str
        :type target_coefficient: float
        :type fraction_of_optimum: float
        :type is_max: bool
        :type pathway_id: str
        :type objective_id: str

        :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
        :rtype: tuple
        """
        #retreive the biomass objective and flux results and set as maxima
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        self.logger.debug('findCreateObjective: '+str(source_reaction))
        source_obj_id = self.findCreateObjective([source_reaction], [source_coefficient], is_max)
        #TODO: use the rpSBML BRSynth annotation parser
        source_flux = None
        try:
            fbc_obj = fbc_plugin.getObjective(source_obj_id)
            #TODO: if this is None need to set it up 
            fbc_obj_annot = fbc_obj.getAnnotation()
            if not fbc_obj_annot:
                raise ValueError
            source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
            self.logger.debug('Already calculated flux for '+str(source_obj_id))
        except (AttributeError, ValueError) as e:
            self.logger.debug('Performing FBA to calculate the source reaction: '+str(source_reaction))
            ### FBA ###
            #self.runFBA(source_reaction, source_coefficient, is_max, pathway_id)
            self._checklibSBML(fbc_plugin.setActiveObjectiveId(source_obj_id),
                    'Setting active objective '+str(source_obj_id))
            self.cobra_model = self._convertToCobra()
            #convert model to cobra since you have made some changes 
            if not self.cobra_model:
                self._writeAnalysisResults(source_obj_id, 0.0, pathway_id)
                return 0.0, False
            cobra_results = self.cobra_model.optimize()
            self.logger.debug('Source reaction optimisation results: '+str(cobra_results.objective_value))
            self._writeAnalysisResults(source_obj_id, cobra_results, pathway_id)
            # cobra_results.objective_value
            fbc_obj = fbc_plugin.getObjective(source_obj_id)
            fbc_obj_annot = fbc_obj.getAnnotation()
            try:
                self._checklibSBML(fbc_obj_annot, 'Retreiving the annotation: '+str(source_obj_id))
            except AttributeError:
                self.logger.error('No annotation available for: '+str(source_obj_id))
                return 0.0, False
            source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
            self.logger.debug('The source flux is: '+str(source_flux))
        #TODO: add another to check if the objective id exists
        self.logger.debug('FBA source flux ('+str(source_reaction)+') is: '+str(source_flux))
        if not objective_id:
            objective_id = 'obj_'+str(target_reaction)+'__restricted_'+str(source_reaction)
        #self.logger.debug('findCreateObjective() for '+str(objective_id))
        objective_id = self.findCreateObjective([target_reaction], [target_coefficient], is_max, objective_id)
        self.logger.debug('Optimising the objective: '+str(objective_id))
        self.logger.debug('Setting upper bound of source reaction ('+str(source_reaction)+'): '+str(source_flux*fraction_of_source))
        self.logger.debug('Setting lower bound of source reaction ('+str(source_reaction)+'): '+str(source_flux*fraction_of_source))
        #### write the restritions the heterologous pathway annotations
        groups = self.model.getPlugin('groups')
        self._checklibSBML(groups, 'Getting groups plugin')
        try:
            rp_pathway = groups.getGroup(pathway_id)
            self._checklibSBML(rp_pathway, 'retreiving the group: '+str(pathway_id))
        except AttributeError:
            self.logger.warning('The group '+str(pathway_id)+' does not exist... creating it')
            self.createPathway(pathway_id)
            rp_pathway = groups.getGroup(pathway_id)
            self._checklibSBML(rp_pathway, 'Getting RP pathway')
        #write the results to the rp_pathway
        self.logger.debug('Adding annotation to '+str(pathway_id)+' for the restriction value: '+str(source_flux*fraction_of_source))
        self.addUpdateBRSynth(rp_pathway, 'fba_obj_'+str(source_reaction)+'_restricted', str(source_flux*fraction_of_source), 'mmol_per_gDW_per_hr', False)
        #ser the restrictions
        old_upper_bound, old_lower_bound = self.setReactionConstraints(source_reaction,
                                                                       source_flux*fraction_of_source,
                                                                       source_flux*fraction_of_source)
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        #convert model to cobra since you have made some changes 
        self.cobra_model = self._convertToCobra()
        if not self.cobra_model:
            self.logger.error('Converting libSBML to CobraPy returned False')
            #although this may not be the greatest idea, set flux to 0.0 when cobrapy error
            self._writeAnalysisResults(objective_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = self.cobra_model.optimize()
        self._writeAnalysisResults(objective_id, cobra_results, pathway_id)
        ##### print the biomass results ######
        #self.logger.debug('Biomass: '+str(cobra_results.fluxes.biomass))
        #self.logger.debug('Target: '+str(cobra_results.fluxes.RP1_sink))
        #reset the bounds to the original values for the target
        old_upper_bound, old_lower_bound = self.setReactionConstraints(source_reaction,
                                                                       old_upper_bound,
                                                                       old_lower_bound)
        self.logger.debug('The objective '+str(objective_id)+' results '+str(cobra_results.objective_value))
        return cobra_results.objective_value, True


    ########################################################################
    ############################### FBA pathway ranking ####################
    ########################################################################

    #1) Number of interventions
    # need to calculate the number of steps that are not native to know the number of interventions

    #2) Maximal growth rate

    #3) Minimum product yeild at maximal growth rate

    #4) Minimum product yeild

    #5) Anaerobic condition

    #6) Number of potentially disruptive products

        #Toxicity?

    #7) Number of accessible metabolites (avoid intermediate accumulation)

    #8) Thermodynamics (MDF)

    #9) The overlap of the same changes --> might not be applicable in our case

    #10) Reduced model

    #11) ECM
