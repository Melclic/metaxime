import tempfile
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


__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = ["Joan Herisson"]
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


class rpSelenzyme(rpSBML):
    """Class to handle calling the Selenzyme service
    """
    def __init__(self,
                 cache_tar_path=os.path.join('metaxime', 'input_cache', 'rpselenzyme_data.tar.xz'),
                 uniprot_aaLength=None,
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
        self.logger = logging.getLogger(__name__)
        self.cache_tar_path = cache_tar_path
        self.uniprot_aaLength = uniprot_aaLength
        self.data_dir = data_dir
        self.pc = pc
        if self.cache_tar_path:
            if not os.path.exists(self.cache_tar_path):
                self.logger.warning('Cannot find the cache tar file: '+str(self.cache_tar_path))
                self.cache_tar_path = None
        if not self.pc or not self.uniprot_aaLength or not self.data_dir:
            self.pc, self.uniprot_aaLength, self.data_dir = self.loadCache(cache_tar_path)
        else:
            self.logger.warning('Either both pc and cache_tar_path are None and the cache is loaded or both are passed')


    #overwrite the del class to remove the self.data_dir if it has been initiated
    def __del__(self):
        self.logger.debug('Destructor called')
        if isinstance(self.data_dir, tempfile.TemporaryDirectory):
            self.logger.debug('Cleaning up: '+str(self.data_dir.name))
            self.data_dir.cleanup()
        

    def loadCache(self, rpselen_cache_tar=os.path.join('metaxime', 'input_cache', 'rpselenzyme_data.tar.xz')):
        """Load the selenzyme cache
        """
        ## global parameter 
        if not os.path.exists(rpselen_cache_tar):
            self.logger.error('The input cache does not seem to exists: '+str(rpselen_cache_tar))
            return None, None, None
        if not os.path.isdir(os.path.join('metaxime', 'cache')):
            os.mkdir(os.path.join('metaxime', 'cache')) 
        pc = None
        uniprot_aaLength = {}
        tmp_output_folder = tempfile.TemporaryDirectory()
        tar = tarfile.open(rpselen_cache_tar, mode='r')
        tar.extractall(path=tmp_output_folder.name)
        tar.close()
        pc = Selenzy.readData(os.path.join(tmp_output_folder.name, 'data'))
        #### load the uniprot_aaLength
        if os.path.exists(os.path.join('metaxime', 'cache', 'uniprot_aaLength.pickle.gz')):
            uniprot_aaLength = pickle.load(gzip.open(os.path.join('metaxime', 'cache', 'uniprot_aaLength.pickle.gz'), 'rb'))
        else: 
            with open(os.path.join(tmp_output_folder.name, 'data', 'sel_len.csv')) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                next(csv_reader)
                for row in csv_reader:
                    uniprot_aaLength[row[0].split('|')[1]] = int(row[1])
            pickle.dump(uniprot_aaLength, gzip.open(os.path.join('metaxime', 'cache', 'uniprot_aaLength.pickle.gz'), 'wb'))
        return pc, uniprot_aaLength, tmp_output_folder


    def singleReactionRule(self,
                           reaction_rule,
                           host_taxonomy_id,
                           num_results=50,
                           direction=0,
                           noMSA=True,
                           fp='RDK',
                           rxntype='smarts'):
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
            success, results = Selenzy.analyse(['-'+rxntype, reaction_rule], num_results, os.path.join(self.data_dir.name, 'data')+'/', tmp_output_folder, 'tmp.csv', 0, host_taxonomy_id, pc=self.pc, NoMSA=noMSA)
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
                            if self.uniprot_aaLength[uniprot]>int(min_aa_length):
                                uniprotID_score_restricted[uniprot] = uniprotID_score[uniprot]
                        except KeyError:
                            self.logger.warning('Cannot find the following UNIPROT '+str(uniprot)+' in uniprot_aaLength')
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
