import tempfile
import time
import glob
import pickle
import tarfile
import gzip
import logging
import json
import sys
import csv
import os

from selenzy import Selenzy
from .rpSBML import rpSBML
from .rpCache import rpCache


__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = ["Joan Herisson"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Melchior du Lac"
__status__ = "Development"



class rpSelenzyme(rpSBML):
    """Class to handle calling the Selenzyme service
    """
    def __init__(self,
                 cache_tar_path=None,
                 uniprot_aa_length=None,
                 data_dir=None,
                 pc=None,
                 model_name=None,
                 document=None,
                 path=None,
                 rpcache=None):
        """Class constructor

        :param pc: The cache of Selenzyme preLoad object


        #WARNING: you have to close the data_dir when you are done!!!!
        """
        super().__init__(model_name, document, path, rpcache)
        #self.logger = logging.getLogger(__name__)
        self.logger = logging.getLogger(os.path.basename(__file__))
        self.dirname = os.path.dirname(os.path.abspath( __file__ ))
        if not cache_tar_path:
            self.cache_tar_path = os.path.join(self.dirname, 'input_cache', 'rpselenzyme_data.tar.xz')
        else:
            self.cache_tar_path = cache_tar_path
        self.logger.debug('cache_tar_path: '+str(self.cache_tar_path))
        self.uniprot_aa_length = uniprot_aa_length
        self.data_dir = data_dir
        self.pc = pc
        if not self.uniprot_aa_length or not self.data_dir or not self.pc:
            if not os.path.exists(self.cache_tar_path):
                self.logger.warning('Cannot find the cache tar file: '+str(self.cache_tar_path))
                self.cache_tar_path = None
            if not self.pc or not self.uniprot_aa_length or not self.data_dir:
                self.pc, self.uniprot_aa_length, self.data_dir = self.loadCache(self.cache_tar_path)
            else:
                self.logger.warning('Either both pc and cache_tar_path are None and the cache is loaded or both are passed')
        else:
            self.logger.warning('Cache has been passed fully')
            self.logger.debug('self.data_dir.name: '+str(self.data_dir.name))


    #overwrite the del class to remove the self.data_dir if it has been initiated
    def __del__(self):
        self.logger.debug('Destructor called')
        if isinstance(self.data_dir, tempfile.TemporaryDirectory):
            self.logger.debug('Cleaning up: '+str(self.data_dir.name))
            self.data_dir.cleanup()


    #############################################################
    #################### STATIC #################################
    #############################################################


    @staticmethod
    def runCollection(rpcollection,
                      host_taxonomy_id,
                      rpcollection_output=None,
                      num_results=50,
                      direction=0,
                      noMSA=True,
                      fp='RDK',
                      rxntype='smarts',
                      min_aa_length=100,
                      pathway_id='rp_pathway',
                      pc=None,
                      uniprot_aa_length=None,
                      cache_path=None,
                      rpcache=None):
        with tempfile.TemporaryDirectory() as tmp_folder:
            tar = tarfile.open(rpcollection, mode='r')
            #get the root member
            root_name = os.path.commonprefix(tar.getnames())
            tar.extractall(path=tmp_folder, members=tar.members)
            tar.close()
            if len(glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')))==0:
                logging.error('Input collection has no models')
                return False
            ##### log ########
            rpselenzyme_log = None
            rpfba_log = None
            if os.path.exists(os.path.join(tmp_folder, root_name, 'log.json')):
                rpselenzyme_log = json.load(open(os.path.join(tmp_folder, root_name, 'log.json')))
            else:
                logging.warning('The log does not seem to exists, creating it...')
                rpselenzyme_log = {}
            if not 'rpselenzyme' in rpselenzyme_log:
                rpselenzyme_log['rpselenzyme'] = {}
            rpselenzyme_log['rpselenzyme'][time.time()] = {'rpcollection': rpcollection,
                                                           'host_taxonomy_id': host_taxonomy_id,
                                                           'rpcollection_output': rpcollection_output,
                                                           'num_results':num_results,
                                                           'direction': direction,
                                                           'noMSA': noMSA,
                                                           'fp': fp,
                                                           'rxntype': rxntype,
                                                           'min_aa_length': min_aa_length,
                                                           'cache_path': cache_path}
            json.dump(rpselenzyme_log, open(os.path.join(tmp_folder, root_name, 'log.json'), 'w'))
            cache_selenzyme = None
            ##### cache ######
            if not rpcache:
                rpcache = rpCache()
                #rpcache.populateCache()
            if not pc or not uniprot_aa_length:
                if not cache_path:
                    cache_selenzyme = rpSelenzyme(cache_tar_path=os.path.join(os.path.dirname(os.path.abspath( __file__ )), 'input_cache', 'rpselenzyme_data.tar.xz'))
                else:
                    cache_selenzyme = rpSelenzyme(cache_tar_path=cache_path)
            for rpsbml_path in glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')):
                file_name = rpsbml_path.split('/')[-1].replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '').replace('_rpsbml', '')
                logging.debug('################### '+str(file_name)+' ##################')
                cache_selenzyme.model_name = file_name
                cache_selenzyme.rpsbml_path = rpsbml_path
                cache_selenzyme.readSBML(rpsbml_path, file_name)
                cache_selenzyme.run(host_taxonomy_id=int(host_taxonomy_id),
                                    pathway_id=pathway_id,
                                    num_results=num_results,
                                    direction=direction,
                                    noMSA=noMSA,
                                    fp=fp,
                                    rxntype=rxntype,
                                    min_aa_length=min_aa_length,
                                    write_results=True)
                cache_selenzyme.writeSBML(path=rpsbml_path)
            if len(glob.glob(os.path.join(tmp_folder, root_name, 'models', '*')))==0:
                logging.error('Output has not produced any models')
                return False
            if isinstance(cache_path, tempfile.TemporaryDirectory):
                cache_path.cleanup()
            #WARNING: we are overwriting the input file
            if rpcollection_output:
                with tarfile.open(rpcollection_output, "w:xz") as tar:
                    tar.add(os.path.join(tmp_folder, root_name), arcname='rpsbml_collection')
            else:
                logging.warning('The output file is: '+str(os.path.join(os.path.dirname(rpcollection), 'output.tar.xz')))
                with tarfile.open(os.path.join(os.path.dirname(rpcollection), 'output.tar.xz'), "w:xz") as tar:
                    tar.add(os.path.join(tmp_folder, root_name), arcname='rpsbml_collection')
        return True 


    #############################################################
    #################### PUBLIC #################################
    #############################################################


    def loadCache(self, cache_tar_path=None):
        """Load the selenzyme cache
        """
        ## global parameter 
        self.logger.debug('cache_tar_path: '+str(cache_tar_path))
        if not cache_tar_path:
            self.logger.error('Cannot have None tar path')
            return None, None, None
        if not os.path.exists(cache_tar_path):
            self.logger.error('The input cache does not seem to exists: '+str(cache_tar_path))
            return None, None, None
        if not os.path.isdir(os.path.join(self.dirname, 'cache')):
            os.mkdir(os.path.join(self.dirname, 'cache')) 
        pc = None
        uniprot_aa_length = {}
        tmp_output_folder = tempfile.TemporaryDirectory()
        tar = tarfile.open(cache_tar_path, mode='r')
        tar.extractall(path=tmp_output_folder.name)
        tar.close()
        pc = Selenzy.readData(os.path.join(tmp_output_folder.name, 'data'))
        #### load the uniprot_aa_length
        if os.path.exists(os.path.join(self.dirname, 'cache', 'uniprot_aa_length.pickle.gz')):
            uniprot_aa_length = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', 'uniprot_aa_length.pickle.gz'), 'rb'))
        else: 
            with open(os.path.join(tmp_output_folder.name, 'data', 'sel_len.csv')) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                next(csv_reader)
                for row in csv_reader:
                    uniprot_aa_length[row[0].split('|')[1]] = int(row[1])
            pickle.dump(uniprot_aa_length, gzip.open(os.path.join(self.dirname, 'cache', 'uniprot_aa_length.pickle.gz'), 'wb'))
        return pc, uniprot_aa_length, tmp_output_folder


    def singleReactionRule(self,
                           reaction_rule,
                           host_taxonomy_id,
                           num_results=50,
                           direction=0,
                           noMSA=True,
                           fp='RDK',
                           rxntype='smarts'): #WARNING: we are actually inputting Reaction SMILES, but it seems Selenzyme made a mistake here
        """Query Selenzyme given a single reaction rule

        :param reaction_rule: The reaction rule
        :param host_taxonomy_id: The taxonomy id associated with the reaction rule
        :param num_results: The number of uniprot ids to return (Default: 50)
        :param direction: Forward (1) to reverse (0) direction (Default: 0)
        :param noMSA: Perform sequence alignement or not (Default: True)
        :param fp: Fingerprint for reactants for quickRSiml (Default: RDK)
        :param rxntype: The type of reaction rule. Valid options: smarts, smiles. (Default: smarts)

        :type reaction_rule: str
        :type host_taxonomy_id: int
        :type num_results: int
        :type direction: int
        :type noMSA: bool
        :type fp: str
        :type rxntype: str

        :rtype: dict
        :return: The UNIPROT id's and its associated score
        """
        uniprotID_score = {}
        score = Selenzy.seqScore()
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            self.logger.debug(os.path.join(self.data_dir.name, 'data'))
            self.logger.debug(os.listdir(os.path.join(self.data_dir.name, 'data')))
            self.logger.debug(tmp_output_folder)
            success, results = Selenzy.analyse(['-'+rxntype, reaction_rule], num_results, os.path.join(self.data_dir.name, 'data'), tmp_output_folder, 'tmp.csv', 0, host_taxonomy_id, pc=self.pc, NoMSA=noMSA)
            self.logger.debug(os.listdir(os.path.join(self.data_dir.name, 'data')))
            self.logger.debug(os.listdir(tmp_output_folder))
            self.logger.debug(os.path.join(tmp_output_folder, 'tmp.csv'))
            self.logger.debug(open(os.path.join(tmp_output_folder, 'tmp.csv'), 'r').read())
            data = Selenzy.updateScore(os.path.join(tmp_output_folder, 'tmp.csv'), score)
            val = json.loads(data.to_json())
            if 'Seq. ID' in val and len(val['Seq. ID'])>0:
                for ix in sorted(val['Seq. ID'], key=lambda z: int(z)):
                    uniprotID_score[val['Seq. ID'][ix]] = val['Score'][ix]
            else:
                raise ValueError
            return uniprotID_score


    def run(self,
            host_taxonomy_id,
            pathway_id='rp_pathway',
            num_results=50,
            direction=0,
            noMSA=True,
            fp='RDK',
            rxntype='smarts',
            min_aa_length=100,
            write_results=False):
        """Return UNIPROT id's associated with each reaction using Selenzyme, from a rpSBML object

        :param host_taxonomy_id: The taxonomy id of the host (Default: 83333, ie. E.Coli)
        :param pathway_id: Group id of the heterologous pathway (Default: rp_pathway)
        :param num_results: Number of UNIPROT id's to return per reaction rule (Default: 50)
        :param direction: Forward (1) to reverse (0) direction (Default: 0)
        :param noMSA: Perform sequence alignement or not (Default: True)
        :param fp: Fingerprint for reactants for quickRSiml (Default: RDK)
        :param rxntype: The type of reaction rule. Valid options: smarts, smiles. (Default: smarts)
        :param min_aa_length: Filter the UNIRPOT proteins and return only whose amino acid lengths are greater than the input value. (Default: 100)

        :type host_taxonomy_id: int
        :type pathway_id: str
        :type num_results: int
        :type direction: int
        :type noMSA: bool
        :type fp: str
        :type rxntype: str
        :type min_aa_length: int

        :rtype: bool
        :return: The success or failure of the function
        """
        #If none are passed, then try to recover from the current model, i.e. they are merged
        if not host_taxonomy_id:
            res = self.model.readTaxonomy()
            if not res:
                self.logger.error('Cannot determine the taxonomy of the input. Please input one')
            else:
                if len(res)>1:
                    self.logger.warning('Multiple returned taxonomy, selecting the first one arbitrarily: '+str(res))
                try:
                    host_taxonomy_id = int(res[0])
                except ValueError:
                    self.logger.error('Cannot convert the following taxonomy id to int: '+str(res[0]))
                    return False
        all_results = {}
        for reac_id in self.getGroupsMembers(pathway_id):
            reac = self.model.getReaction(reac_id)
            brs_reac = self.readBRSYNTHAnnotation(reac.getAnnotation())
            all_results[reac_id] = {}
            if brs_reac['smiles']:
                try:
                    uniprotID_score = self.singleReactionRule(brs_reac['smiles'], host_taxonomy_id, num_results, direction, noMSA, fp, rxntype)
                    uniprotID_score_restricted = {}
                    for uniprot in uniprotID_score:
                        try:
                            if self.uniprot_aa_length[uniprot]>int(min_aa_length):
                                uniprotID_score_restricted[uniprot] = uniprotID_score[uniprot]
                        except KeyError:
                            self.logger.warning('Cannot find the following UNIPROT '+str(uniprot)+' in uniprot_aa_length')
                    xref = {'uniprot': [i for i in uniprotID_score_restricted]}
                    if write_results:
                        self.addUpdateMIRIAM(reac, 'reaction', xref)
                        self.addUpdateBRSynth(reac, 'selenzyme', uniprotID_score_restricted, None, False, True, True)
                    all_results[reac_id] = uniprotID_score_restricted
                except ValueError:
                    self.logger.warning('Problem with retreiving the selenzyme information for model '+str(self.model.getId()))
                    return False
            else:
                self.logger.warning('Cannot retreive the reaction rule of model '+str(self.model.getId()))
                return False
        return all_results
