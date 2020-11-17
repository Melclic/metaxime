import libsbml
import sys
import os
import json
import tempfile
import glob
import os
import tarfile
import logging
import cobra
import time
from cobra.flux_analysis import pfba
from multiprocessing import Pool

from shared_memory_dict import SharedMemoryDict

from .rpMerge import rpMerge
from .rpCache import rpCache


__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = []
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Melchior du Lac"
__status__ = "Development"

logging.root.setLevel(logging.NOTSET)

logging.basicConfig(
    #level=logging.DEBUG,
    level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


###########################################################
######################### main ############################
###########################################################


#TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes
class rpFBA(rpMerge):
    """Class to simulate an rpsbml object using different FBA types and objective functions
    """
    def __init__(self,
                 sbml_path=None,
                 is_gem_sbml=False,
                 model_name=None,
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
                         path=sbml_path,
                         rpcache=rpcache)
        self.logger = logging.getLogger(__name__)
        self.logger.debug('Started instance of rpFBA')
        self.species_source_target = None
        self.reactions_source_target = None
        self.sbml_path = sbml_path
        self.is_gem_sbml = is_gem_sbml
        self.model_name = model_name
        self.rpcache = rpcache
        self.cobra_model = None


    ##########################################################
    ################## Static Functions ######################
    ##########################################################


    @staticmethod
    def multiFBA(gem_sbml_path,
                 rpsbml_path,
                 shared_rpcache_name,
                 model_name=None,
                 del_sp_pro=False,
                 del_sp_react=True,
                 upper_flux_bound=999999.0,
                 lower_flux_bound=0.0,
                 compartment_id='MNXC3',
                 pathway_id='rp_pathway',
                 created_reaction_pathway_id='merge_created_reactions',
                 sim_type='fraction',
                 keep_merged=False,
                 source_reaction='biomass',
                 source_coefficient=1.0,
                 target_reaction='RP1_sink',
                 target_coefficient=1.0,
                 fraction_of=0.75,
                 is_max=True,
                 objective_id='obj_fraction'):
        logging.debug('################ '+str(model_name)+' ###############')
        rpfba = rpFBA(sbml_path=gem_sbml_path,
                      is_gem_sbml=True,
                      model_name=model_name,
                      rpcache=None)
        """Function to run as process
        """
        if shared_rpcache_name:
            rpcache_shared_dict = SharedMemoryDict(name=shared_rpcache_name, size=None) #not sure why I need to pass size although not used
            rpfba.setFromDict(rpcache_shared_dict)
            '''
            logging.debug('rpcache_shared_dict[deprecatedCID_cid]: '+str(hex(id(rpcache_shared_dict['deprecatedCID_cid']))))
            logging.debug('rpcache_shared_dict[deprecatedCID_cid][MNXM95986]: '+str(rpcache_shared_dict['deprecatedCID_cid']['MNXM95986']))
            logging.debug('rpcache_shared_dict[deprecatedCID_cid][MNXM95986]: '+str(hex(id(rpcache_shared_dict['deprecatedCID_cid']['MNXM95986']))))
            logging.debug('rpfba.deprecatedCID_cid: '+str(hex(id(rpfba.deprecatedCID_cid))))
            logging.debug('rpfba.deprecatedCID_cid[MNXM95986]: '+str(rpfba.deprecatedCID_cid['MNXM95986']))
            logging.debug('rpfba.deprecatedCID_cid[MNXM95986]: '+str(hex(id(rpfba.deprecatedCID_cid['MNXM95986']))))
            '''
        else:
            logging.warning('No shared memory: rpcache_shared_dict')
        status = rpfba.mergeModels(input_model=rpsbml_path,
                                   del_sp_pro=del_sp_pro,
                                   del_sp_react=del_sp_react,
                                   upper_flux_bound=upper_flux_bound,
                                   lower_flux_bound=lower_flux_bound,
                                   compartment_id=compartment_id,
                                   pathway_id=pathway_id,
                                   created_reaction_pathway_id=created_reaction_pathway_id)
        if not status:
            logging.error('Problem merging the models: '+str(file_name))
            return False
        ####### fraction of reaction ######
        if sim_type=='fraction':
            status = rpfba.runFractionReaction(source_reaction, source_coefficient, target_reaction, target_coefficient, fraction_of, is_max, pathway_id, objective_id, True)
        ####### FBA ########
        elif sim_type=='fba':
            status = rpfba.runFBA(target_reaction, target_coefficient, is_max, pathway_id, objective_id, True)
        #TODO: add multi objective
        ####### pFBA #######
        elif sim_type=='pfba':
            status = rpfba.runParsimoniousFBA(target_reaction, target_coefficient, fraction_of, is_max, pathway_id, objective_id, True)
        if not status:
            logging.error('Problem running the FBA')
            return False
        #cleanup the shared memory
        if shared_rpcache_name:
            rpcache_shared_dict.cleanup()
        #overwrite the rpSBML model
        if keep_merged:
            return rpfba.writeSBML(path=rpsbml_path)
        else:
            return rpfba.writeSBML(path=rpsbml_path, skinny_rpsbml_path=rpsbml_path)


    @staticmethod
    def runCollection(rpcollection,
                      gem_sbml_path,
                      rpcollection_output=None,
                      target_reaction='biomass',
                      sim_type='fba',
                      source_reaction=None,
                      source_coefficient=1.0,
                      target_coefficient=1.0,
                      num_workers=10,
                      keep_merged=False,
                      fraction_of=0.75,
                      objective_id='obj_fba',
                      is_max=True,
                      del_sp_pro=False,
                      del_sp_react=True,
                      upper_flux_bound=999999.0,
                      lower_flux_bound=0.0,
                      compartment_id='MNXC3',
                      pathway_id='rp_pathway',
                      created_reaction_pathway_id='merge_created_reactions',
                      rpcache=None):
        """Parse a the collection of rpSBML files and its extra information

        Note: That this function does nothing. It is written to be extended from the other inherited classes
        TODO: Not sure of the rpcache passing, only applicable when you are extending parsing many different

        """
        ######## main #######
        if not os.path.exists(gem_sbml_path):
            logging.error('The input GEM file does not seem to exists: '+str(gem_sbml_path))
            return False
        if not sim_type in ['fraction', 'fba', 'pfba']:
            logging.error('The input FBA sim_types is not supported: '+str(sim_type))
            return False
        if num_workers<1 or num_workers>50:
            logging.error('Cannot have 0 or more than 50 workers')
            return False
        with tempfile.TemporaryDirectory() as tmp_folder:
            tar = tarfile.open(rpcollection, mode='r')
            #get the root member
            root_name = os.path.commonprefix(tar.getnames())
            logging.debug('root_name: '+str(root_name))
            tar.extractall(path=tmp_folder, members=tar.members)
            tar.close()
            #check and create the log
            rpfba_log = None
            if os.path.exists(os.path.join(tmp_folder, root_name, 'log.json')):
                rpfba_log = json.load(open(os.path.join(tmp_folder, root_name, 'log.json')))
            else:
                logging.warning('The log does not seem to exists, creating it...')
                rpfba_log = {}
            if not 'rpfba' in rpfba_log:
                rpfba_log['rpfba'] = {}
            rpfba_log['rpfba'][time.time()] = {'gem_sbml_path': gem_sbml_path,
                                               'rpcollection_output': rpcollection_output,
                                               'target_reaction': target_reaction,
                                               'sim_type': sim_type,
                                               'source_reaction': source_reaction,
                                               'source_coefficient': source_coefficient,
                                               'target_coefficient': target_coefficient,
                                               'num_workers': num_workers,
                                               'keep_merged': keep_merged,
                                               'fraction_of': fraction_of,
                                               'objective_id': objective_id,
                                               'is_max': is_max,
                                               'del_sp_pro': del_sp_pro,
                                               'del_sp_react': del_sp_react,
                                               'upper_flux_bound': upper_flux_bound,
                                               'lower_flux_bound': lower_flux_bound,
                                               'compartment_id': compartment_id,
                                               'created_reaction_pathway_id': created_reaction_pathway_id,
                                               'pathway_id': pathway_id}
            json.dump(rpfba_log, open(os.path.join(tmp_folder, root_name, 'log.json'), 'w'))
            if len(glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')))==0:
                logging.error('Input collection has no models')
                return False
            #load the cache
            if not rpcache:
                #TODO: check what part of the cache is used here and only populate that part
                rpcache = rpCache()
                rpcache.populateCache()
            #if only a single worker then we do not use processify -- easier debugging
            if num_workers==1:
                for rpsbml_path in glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')):
                    file_name = rpsbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '')
                    logging.debug('################ '+str(file_name)+' ###############')
                    rpfba = rpFBA(sbml_path=gem_sbml_path,
                                  is_gem_sbml=True,
                                  model_name=file_name,
                                  rpcache=rpcache)
                    status = rpfba.mergeModels(input_model=rpsbml_path,
                                               del_sp_pro=del_sp_pro,
                                               del_sp_react=del_sp_react,
                                               upper_flux_bound=upper_flux_bound,
                                               lower_flux_bound=lower_flux_bound,
                                               compartment_id=compartment_id,
                                               pathway_id=pathway_id,
                                               created_reaction_pathway_id=created_reaction_pathway_id)
                    #debug
                    #logging.debug('write sbml: {0}'.format(rpfba.writeSBML(path='/home/mdulac/Downloads/test.sbml')))
                    if not status:
                        logging.error('Problem merging the models: '+str(file_name))
                        continue
                    ####### fraction of reaction ######
                    if sim_type=='fraction':
                        rpfba.runFractionReaction(source_reaction, source_coefficient, target_reaction, target_coefficient, fraction_of, is_max, pathway_id, objective_id, True)
                    ####### FBA ########
                    elif sim_type=='fba':
                        rpfba.runFBA(target_reaction, target_coefficient, is_max, pathway_id, objective_id, True)
                    ####### pFBA #######
                    elif sim_type=='pfba':
                        rpfba.runParsimoniousFBA(target_reaction, target_coefficient, fraction_of, is_max, pathway_id, objective_id, True)
                    #overwrite the rpSBML model
                    if keep_merged:
                        rpfba.writeSBML(path=rpsbml_path)
                    else:
                        rpfba.writeSBML(path=rpsbml_path, skinny_rpsbml_path=rpsbml_path)
            else:
                logging.debug('Setting up the shared memory')
                rpcache_shared_dict = SharedMemoryDict(name='rpcache', size=rpcache.getSizeCache())
                rpcache_shared_dict['cid_strc'] = rpcache.cid_strc
                rpcache_shared_dict['deprecatedCID_cid'] = rpcache.deprecatedCID_cid
                rpcache_shared_dict['deprecatedRID_rid'] = rpcache.deprecatedRID_rid
                rpcache_shared_dict['cid_xref'] = rpcache.cid_xref
                rpcache_shared_dict['comp_xref'] = rpcache.comp_xref
                rpcache_shared_dict['rr_full_reactions'] = rpcache.rr_full_reactions
                rpcache_shared_dict['chebi_cid'] = rpcache.chebi_cid
                rpcache_shared_dict['inchikey_cid'] = rpcache.inchikey_cid
                rpcache_shared_dict['cid_name'] = rpcache.cid_name
                #lock the shared memory to read only
                result_objs = []
                with Pool(processes=num_workers) as pool:
                    for rpsbml_path in glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')):
                        file_name = rpsbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '')
                        result = pool.apply_async(rpFBA.multiFBA, args=(gem_sbml_path,
                                                                        rpsbml_path,
                                                                        'rpcache',
                                                                        file_name,
                                                                        del_sp_pro,
                                                                        del_sp_react,
                                                                        upper_flux_bound,
                                                                        lower_flux_bound,
                                                                        compartment_id,
                                                                        pathway_id,
                                                                        created_reaction_pathway_id,
                                                                        sim_type,
                                                                        keep_merged,
                                                                        source_reaction,
                                                                        source_coefficient,
                                                                        target_reaction,
                                                                        target_coefficient,
                                                                        fraction_of,
                                                                        is_max,
                                                                        objective_id,))
                        result_objs.append(result)
                        logging.debug('result: '+str(result))
                    results = [result.get() for result in result_objs]
                    logging.debug('results: '+str(results))
                #clean up the 
                rpcache_shared_dict.cleanup()
            if len(glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')))==0:
                logging.error('Output has not produced any models')
                return False
            #WARNING: we are overwriting the input file
            if rpcollection_output:
                with tarfile.open(rpcollection_output, "w:xz") as tar:
                    tar.add(os.path.join(tmp_folder, root_name), arcname='rpsbml_collection')
            else:
                logging.warning('The output file is: '+str(os.path.join(os.path.dirname(rpcol), 'output.tar.xz')))
                with tarfile.open(os.path.join(os.path.dirname(rpcol), 'output.tar.xz'), "w:xz") as tar:
                    tar.add(os.path.join(tmp_folder, root_name), arcname='rpsbml_collection')
        return True


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
                self.writeSBML(os.path.join(tmp_output_folder, 'tmp.sbml'))
                #logging.info(glob.glob(tmp_output_folder+'/*'))
                #logging.info(cobra.io.validate_sbml_model(glob.glob(tmp_output_folder+'/*')[0]))
                cobra_model = cobra.io.read_sbml_model(os.path.join(tmp_output_folder, 'tmp.sbml'), use_fbc_package=True)
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
        if isinstance(cobra_results, float):
            obj_value = cobra_results 
        else:
            obj_value = cobra_results.objective_value
        self.logger.debug('Set '+str(pathway_id)+' with '+str('fba_'+str(objective_id))+' to '+str(obj_value))
        self.addUpdateBRSynth(rp_pathway, 'fba_'+str(objective_id), str(obj_value), 'mmol_per_gDW_per_hr', False)
        #get the objective
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC plugin')
        obj = fbc_plugin.getObjective(objective_id)
        self._checklibSBML(obj, 'Getting objective '+str(objective_id))
        self.addUpdateBRSynth(obj, 'flux_value', str(obj_value), 'mmol_per_gDW_per_hr', False)
        self.logger.debug('Set the objective '+str(objective_id)+' a flux_value of '+str(obj_value))
        #these are the reaction ones
        if not isinstance(cobra_results, float):
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
                sol_flux_dict = cobra_results.to_frame().to_dict()
                reaction = self.model.getReaction(member.getIdRef())
                self._checklibSBML(reaction, 'Fetching the reaction: '+str(member.getIdRef()))
                try:
                    member_flux = sol_flux_dict['fluxes'][member.getIdRef()]
                except KeyError:
                    self.logger.warning('Cannot find the following reaction in cobrapy fluxes solution: '+str(member.getIdRef()))
                    continue
                self.logger.debug('Set the reaction '+str(member.getIdRef())+' a '+str('fba_'+str(objective_id))+' of '+str(member_flux))
                self.addUpdateBRSynth(reaction, 'fba_'+str(objective_id), str(member_flux), 'mmol_per_gDW_per_hr', False)


    ##################################################################
    ######################## PUBLIC FUNCTIONS ########################
    ##################################################################


    def writeSBML(self, path, skinny_rpsbml_path=None):
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
        if skinny_rpsbml_path:
            if not self.path:
                self.logger.error('Must have path initiated to return the non-merged version of the model')
                return False
            #open the original file once again
            skinny_rpmerge = rpMerge(is_gem_sbml=False,
                                     path=skinny_rpsbml_path,
                                     rpcache=self.rpcache)
            #TODO: need to add the newly created reaction from the merge to the skinny rpSBML version
            #overwrite the original BRSynth annotations with the new ones
            self.logger.debug('asDict: '+str(self.asDict()))
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
            if write_results:
                self._writeAnalysisResults(objective_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = self.cobra_model.optimize()
        if write_results:
            self._writeAnalysisResults(objective_id, cobra_results, pathway_id)
        return 0.0, True


    #TODO
    #def runMultiObjectiveParsimonou


    def runFBA(self, reaction_id, coefficient=1.0, is_max=True, pathway_id='rp_pathway', objective_id=None, write_results=False):
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
        self.logger.debug('objective_id: '+str(objective_id))
        objective_id = self.findCreateObjective([reaction_id], [coefficient], is_max, objective_id)
        self.logger.debug('objective_id: '+str(objective_id))
        #run the FBA
        self._checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                'Setting active objective '+str(objective_id))
        if not self.cobra_model:
            self.cobra_model = self._convertToCobra()
        if not self.cobra_model:
            self.logger.error('Cannot convert the model to cobra')
            if write_results:
                self._writeAnalysisResults(objective_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = self.cobra_model.optimize()
        if write_results:
            self._writeAnalysisResults(objective_id, cobra_results, pathway_id)
        return cobra_results.objective_value, True


    def runParsimoniousFBA(self, reaction_id, coefficient=1.0, fraction_of_optimum=0.95, is_max=True, pathway_id='rp_pathway', objective_id=None, write_results=False):
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
        if write_results:
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
                            objective_id='obj_fraction',
                            write_results=False):
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
            if write_results:
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
            if write_results:
                self._writeAnalysisResults(objective_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = self.cobra_model.optimize()
        if write_results:
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
