#!/usr/bin/env python3

import libsbml
import sys
import logging
from scipy import stats
import numpy as np
import json
import glob
import tarfile
import tempfile

from .rpSBML import rpSBML

__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = [""]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Melchior du Lac"
__status__ = "Development"


self.logger.basicConfig(
    #level=self.logger.DEBUG,
    level=self.logger.WARNING,
    #level=self.logger.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


class rpGlobalScore(rpSBML):
    """Class combining all the different characteristics of a pathway and calculate the global score
    """
    def __init__(self
                 model_name=None,
                 document=None,
                 path=None,
                 rpcache=None):
        super().__init__(model_name, document, path, rpcache)
        self.logger = logging.getLogger(__name__)
        self.logger.info('Starting instance of rpGlobalScore')

    '''
    ## Normalise by sigmoidal function
    # NOT USED
    #
    def nonlin(x,deriv=False):
        if(deriv==True):
            return x*(1-x)
        return 1/(1+np.exp(-x))
    '''


    ##################################################################################
    ############################### STATIC ###########################################
    ##################################################################################


    @staticmethod
    def runCollection(rpcollection,
                      rpcollection_output=None,
                      rpcache=None,
                      weight_rp_steps=0.10002239003499142,
                      weight_rule_score=0.13346271414277305,
                      weight_fba=0.6348436269211155,
                      weight_thermo=0.13167126890112002,
                      max_rp_steps=15, #TODO: add this as a limit in RP2
                      thermo_ceil=5000.0,
                      thermo_floor=-5000.0,
                      fba_ceil=5.0,
                      fba_floor=0.0,
                      pathway_id='rp_pathway',
                      objective_id='obj_fraction',
                      thermo_id='dfG_prime_m'):
        with tempfile.TemporaryDirectory() as tmp_folder:
            tar = tarfile.open(rpcollection, mode='r')
            #get the root member
            root_name = os.path.commonprefix(tar.getnames())
            logging.debug('root_name: '+str(root_name))
            tar.extractall(path=tmp_folder, members=tar.members)
            tar.close()  
            if len(glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')))==0:
                logging.error('Input collection has no models')
                return False
            ####### log #########
            rpglobalscore_log = None
            rpfba_log = None
            if os.path.exists(os.path.join(tmp_folder, root_name, 'log.json')):
                rpglobalscore_log = json.load(open(os.path.join(tmp_folder, root_name, 'log.json')))
            else:
                logging.warning('The log does not seem to exists, creating it...')
                rpglobalscore_log = {}
            if not 'rpglobalscore' in rpglobalscore_log:
                rpglobalscore_log['rpglobalscore'] = {}
            rpglobalscore_log['rpglobalscore'][time.time()] = {'rpcollection': rpcollection,
                                                               'rpcollection_output': rpcollection_output,
                                                               'weight_rp_steps': weight_rp_steps,
                                                               'weight_rule_score': weight_rule_score,
                                                               'weight_fba': weight_fba,
                                                               'weight_thermo': weight_thermo,
                                                               'max_rp_steps': max_rp_steps,
                                                               'thermo_ceil': thermo_ceil,
                                                               'thermo_floor': thermo_floor,
                                                               'fba_ceil': fba_ceil,
                                                               'fba_floor': fba_floor,
                                                               'pathway_id': pathway_id,
                                                               'objective_id': objective_id,
                                                               'thermo_id': thermo_id}
            json.dump(rpglobalscore_log, open(os.path.join(tmp_output_folder, root_name, 'log.json'), 'w'))
            for rpsbml_path in glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')):
                file_name = rpsbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '')
                rpglobalscore = rpGlobalScore(model_name=file_name, path=rpsbml_path, rpcache=rpcache)
                rpglobalscore.calculateGlobalScore(weight_rp_steps,
                                                   weight_rule_score,
                                                   weight_fba,
                                                   weight_thermo,
                                                   max_rp_steps2
                                                   thermo_ceil,
                                                   thermo_floor,
                                                   fba_ceil,
                                                   fba_floor,
                                                   pathway_id,
                                                   objective_id,
                                                   thermo_id,
                                                   True)
                rpglobalscore.writeSBML(path=rpsbml_path)
        if len(glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')))==0:
            logging.error('Output has not produced any models')
            return False
        #WARNING: we are overwriting the input file
        if rpcollection_output:
            with tarfile.open(rpcollection_output, "w:xz") as tar:
                tar.add(os.path.join(tmp_folder, root_name), arcname='rpsbml_collection')
        else:
            logging.warning('The output file is: '+str(os.path.join(os.path.dirname(rpcollection), 'output.tar.xz')))
            with tarfile.open(os.path.join(os.path.dirname(rpcollection), 'output.tar.xz'), "w:xz") as tar:
                tar.add(os.path.join(tmp_folder, root_name), arcname='rpsbml_collection')
        return True 


    ##################################################################################
    ############################### PUBLIC ###########################################
    ##################################################################################


    # NOTE: all the scores are normalised by their maximal and minimal, and normalised to be higher is better
    # Higher is better
    #TODO: try to standardize the values instead of normalisation.... Advantage: not bounded
    def calculateGlobalScore(self,
                             weight_rp_steps=0.10002239003499142,
                             weight_rule_score=0.13346271414277305,
                             weight_fba=0.6348436269211155,
                             weight_thermo=0.13167126890112002,
                             max_rp_steps=15, #TODO: add this as a limit in RP2
                             thermo_ceil=5000.0,
                             thermo_floor=-5000.0,
                             fba_ceil=5.0,
                             fba_floor=0.0,
                             pathway_id='rp_pathway',
                             objective_id='obj_fraction',
                             thermo_id='dfG_prime_m',
                             write_results=False):
        """ Extract the reaction SMILES from an SBML, query rule_score and write the results back to the SBML

        :param weight_rp_steps: The weight associated with the number of steps
        :param weight_rule_score: The weight of the mean rule score
        :param weight_fba: The weight of the target flux
        :param weight_thermo: The weight of the thermodynamics
        :param max_rp_steps: The maximum number of steps
        :param thermo_ceil: The maximum pathway thermodynamics value
        :param thermo_floor: The minimum thermodynamics value
        :param fba_ceil: The maximum target fba
        :param fba_floor: The minimum target fba
        :param pathway_id: The id of the heterologous pathway
        :param objective_id: The id if the fba objective id
        :param thermo_id: The is of the thermodynamics

        :type weight_rp_steps: float
        :type weight_rule_score: float
        :type weight_fba: float
        :type weight_thermo: float
        :type max_rp_steps: int
        :type thermo_ceil: float
        :type thermo_floor: float
        :type fba_ceil: float
        :type fba_floor: float
        :type pathway_id: str
        :type objective_id: str
        :type thermo_id: str

        :rtype: float
        :return: The global score
        """
        rpsbml_dict = self.asdict(pathway_id)
        path_norm = {}
        ########################################### REACTIONS #################################################
        #WARNING: we do this because the list gets updated
        self.logger.debug('thermo_ceil: '+str(thermo_ceil))
        self.logger.debug('thermo_floor: '+str(thermo_floor))
        self.logger.debug('fba_ceil: '+str(fba_ceil))
        self.logger.debug('fba_floor: '+str(fba_floor))
        list_reac_id = list(rpsbml_dict['reactions'].keys())
        for reac_id in list_reac_id:
            list_bd_id = list(rpsbml_dict['reactions'][reac_id]['brsynth'].keys())
            for bd_id in list_bd_id:
                ####### Thermo ############
                #lower is better -> -1.0 to have highest better
                #WARNING: we will only take the dfG_prime_m value
                if bd_id[:4]=='dfG_':
                    if bd_id not in path_norm:
                        path_norm[bd_id] = []
                    try:
                        if thermo_ceil>=rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']>=thermo_floor:
                            #min-max feature scaling
                            norm_thermo = (rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']-thermo_floor)/(thermo_ceil-thermo_floor)
                            norm_thermo = 1.0-norm_thermo
                        elif rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']<thermo_floor:
                            norm_thermo = 1.0
                        elif rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']>thermo_ceil:
                            norm_thermo = 0.0
                    except (KeyError, TypeError) as e:
                        self.logger.warning('Cannot find the thermo: '+str(bd_id)+' for the reaction: '+str(reac_id))
                        norm_thermo = 1.0
                    rpsbml_dict['reactions'][reac_id]['brsynth']['norm_'+bd_id] = {}
                    rpsbml_dict['reactions'][reac_id]['brsynth']['norm_'+bd_id]['value'] = norm_thermo
                    self.logger.debug(str(bd_id)+': '+str(rpsbml_dict['reactions'][reac_id]['brsynth']['norm_'+bd_id]['value'])+' ('+str(norm_thermo)+')')
                    path_norm[bd_id].append(norm_thermo)
                ####### FBA ##############
                #higher is better
                #return all the FBA values
                #------- reactions ----------
                elif bd_id[:4]=='fba_':
                    try:
                        norm_fba = 0.0
                        if fba_ceil>=rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']>=fba_floor:
                            #min-max feature scaling
                            norm_fba = (rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']-fba_floor)/(fba_ceil-fba_floor)
                        elif rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']<=fba_floor:
                            norm_fba = 0.0
                        elif rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']>fba_ceil:
                            norm_fba = 1.0
                        rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value'] = norm_fba
                    except (KeyError, TypeError) as e:
                        norm_fba = 0.0
                        self.logger.warning('Cannot find the objective: '+str(bd_id)+' for the reaction: '+str(reac_id))
                    rpsbml_dict['reactions'][reac_id]['brsynth']['norm_'+bd_id] = {}
                    rpsbml_dict['reactions'][reac_id]['brsynth']['norm_'+bd_id]['value'] = norm_fba
                elif bd_id=='rule_score':
                    if bd_id not in path_norm:
                        path_norm[bd_id] = []
                    #rule score higher is better
                    path_norm[bd_id].append(rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value'])
                else:
                    self.logger.debug('Not normalising: '+str(bd_id))
        ########################################### PATHWAY ###################################################
        ############### FBA ################
        #higher is better
        list_path_id = list(rpsbml_dict['pathway']['brsynth'].keys())
        for bd_id in list_path_id:
            if bd_id[:4]=='fba_':
                norm_fba = 0.0
                if fba_ceil>=rpsbml_dict['pathway']['brsynth'][bd_id]['value']>=fba_floor:
                    #min-max feature scaling
                    norm_fba = (rpsbml_dict['pathway']['brsynth'][bd_id]['value']-fba_floor)/(fba_ceil-fba_floor)
                elif rpsbml_dict['pathway']['brsynth'][bd_id]['value']<=fba_floor:
                    norm_fba = 0.0
                elif rpsbml_dict['pathway']['brsynth'][bd_id]['value']>fba_ceil:
                    norm_fba = 1.0
                else:
                    self.logger.warning('This flux event should never happen: '+str(rpsbml_dict['pathway']['brsynth'][bd_id]['value']))
                rpsbml_dict['pathway']['brsynth']['norm_'+bd_id] = {}
                rpsbml_dict['pathway']['brsynth']['norm_'+bd_id]['value'] = norm_fba
        ############# thermo ################
        for bd_id in path_norm:
            if bd_id[:4]=='dfG_':
                rpsbml_dict['pathway']['brsynth']['norm_'+bd_id] = {}
                rpsbml_dict['pathway']['brsynth']['var_'+bd_id] = {}
                #here add weights based on std
                self.logger.debug(str(bd_id)+': '+str(path_norm[bd_id]))
                rpsbml_dict['pathway']['brsynth']['norm_'+bd_id]['value'] = np.average([np.average(path_norm[bd_id]), 1.0-np.std(path_norm[bd_id])], weights=[0.5, 0.5])
                #the score is higher is better - (-1 since we want lower variability)
                #rpsbml_dict['pathway']['brsynth']['var_'+bd_id]['value'] = 1.0-np.var(path_norm[bd_id])
        ############# rule score ############
        #higher is better
        if not 'rule_score' in path_norm:
            self.logger.warning('Cannot detect rule_score: '+str(path_norm))
            rpsbml_dict['pathway']['brsynth']['norm_'+bd_id] = {}
            rpsbml_dict['pathway']['brsynth']['norm_'+bd_id]['value'] = 0.0
        else:
            rpsbml_dict['pathway']['brsynth']['norm_'+bd_id] = {}
            rpsbml_dict['pathway']['brsynth']['norm_'+bd_id]['value'] = np.average(path_norm[bd_id])
        ##### length of pathway ####
        #lower is better -> -1.0 to reverse it
        norm_steps = 0.0
        if len(rpsbml_dict['reactions'])>max_rp_steps:
            self.logger.warning('There are more steps than specified')
            norm_steps = 0.0
        else:
            try:
                norm_steps = (float(len(rpsbml_dict['reactions']))-1.0)/(float(max_rp_steps)-1.0)
                norm_steps = 1.0-norm_steps
            except ZeroDivisionError:
                norm_steps = 0.0
        ################################### GLOBAL ######################################
        ##### global score #########
        try:
            rpsbml_dict['pathway']['brsynth']['norm_steps'] = {}
            rpsbml_dict['pathway']['brsynth']['norm_steps']['value'] = norm_steps
            self.logger.debug('Using the following values for the global score:')
            self.logger.debug('Rule Score: '+str(rpsbml_dict['pathway']['brsynth']['norm_rule_score']['value']))
            self.logger.debug('Thermo: '+str(rpsbml_dict['pathway']['brsynth']['norm_'+str(thermo_id)]['value']))
            self.logger.debug('Steps: '+str(rpsbml_dict['pathway']['brsynth']['norm_steps']['value']))
            self.logger.debug('FBA ('+str('norm_fba_'+str(objective_id))+'): '+str(rpsbml_dict['pathway']['brsynth']['norm_fba_'+str(objective_id)]['value']))
            globalScore = np.average([rpsbml_dict['pathway']['brsynth']['norm_rule_score']['value'],
                                      rpsbml_dict['pathway']['brsynth']['norm_'+str(thermo_id)]['value'],
                                      rpsbml_dict['pathway']['brsynth']['norm_steps']['value'],
                                      rpsbml_dict['pathway']['brsynth']['norm_fba_'+str(objective_id)]['value']],
                                      weights=[weight_rule_score, weight_thermo, weight_rp_steps, weight_fba])
            '''
            globalScore = (rpsbml_dict['pathway']['brsynth']['norm_rule_score']['value']*weight_rule_score+
                           rpsbml_dict['pathway']['brsynth']['norm_'+str(thermo_id)]['value']*weight_thermo+
                           rpsbml_dict['pathway']['brsynth']['norm_steps']['value']*weight_rp_steps+
                           rpsbml_dict['pathway']['brsynth']['norm_fba_'+str(objective_id)]['value']*weight_fba
                           )/sum([weight_rule_score, weight_thermo, weight_rp_steps, weight_fba])
            '''
        except ZeroDivisionError:
            globalScore = 0.0
        except KeyError as e:
            #self.logger.error(rpsbml_dict['pathway']['brsynth'].keys())
            self.logger.error('KeyError for :'+str(e))
            if 'dfG' in e:
                self.logger.error('Have you ran the thermodynamics?')
            elif 'fba' in e:
                self.logger.error('Have you run the FBA on the heterologous pathways?')
            globalScore = 0.0
        rpsbml_dict['pathway']['brsynth']['global_score'] = {}
        rpsbml_dict['pathway']['brsynth']['global_score']['value'] = globalScore
        if write_results:
            self.updateBRSynthPathway(rpsbml_dict, pathway_id)
        return globalScore
