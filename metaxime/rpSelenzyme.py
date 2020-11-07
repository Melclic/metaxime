import tempfile
import pickle
import gzip
import logging
import json
import sys
import csv
import Selenzy
import os


class rpSelenzyme(rpSBML):
    """Class to handle calling the Selenzyme service
    """
    def __init__(self,
                 model_name=None,
                 document=None,
                 path=None,
                 rpcache=None,
                 pc=None,
                 uniprot_aaLenght=None,
                 cache_tar_path='input_cache/rpselenzyme_data.tar.xz'):
        """Class constructor

        :param pc: The cache of Selenzyme preLoad object
        """
        super().__init__(model_name, document, path, rpcache)
        self.logger = logging.getLogger(__name__)
        self.uniprot_aaLenght = None
        self.pc = None
        if (pc and not uniprot_aaLenght) or (not pc and uniprot_aaLenght):
            #self.logger.error('Either both pc and uniprot_aaLenght are None and the cache is loaded or both are passed')
            raise ValueError('Either both pc and uniprot_aaLenght are None and the cache is loaded or both are passed')
        elif uniprot_aaLenght and pc:
            self.uniprot_aaLenght = uniprot_aaLenght
            self.pc = pc
        else:
            if cache_tar_path:
                if os.path.exists(cache_tar_path):
                    self._loadCache()
                else:
                    #self.logger.error('The input cache file is invalid: '+str(cache_tar_path))
                    raise ValueError('The input cache file is invalid: '+str(cache_tar_path))
            else:
                #self.logger.error('Need to specify either the cache tar path, or the Selenzyme pc and uniprot_aaLength')
                raise ValueError('Need to specify either the cache tar path, or the Selenzyme pc and uniprot_aaLength')


    def _loadCache(self, rpselen_cache_tar='input_cache/rpselenzyme_data.tar.xz'):
        """Load the selenzyme cache
        """
        ## global parameter 
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            tar = tarfile.open(rpselen_cache_tar, mode='r')
            tar.extractall(path=tmp_output_folder)
            tar.close()
            self.pc = Selenzy.readData(tmp_output_folder)
            self.uniprot_aaLenght = pickle.load(gzip.open(os.path.join(tmp_output_folder, 'uniprot_aaLenght.pickle.gz'), 'rb'))
            '''This is how the above parameter was generated
            with open(os.path.exists(self.datadir, 'sel_len.csv')) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                next(csv_reader)
                for row in csv_reader:
                    self.uniprot_aaLenght[row[0].split('|')[1]] = int(row[1])
            '''


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
            success, results = Selenzy.analyse(['-'+rxntype, reaction_rule], num_results, self.datadir, tmp_output_folder, 'tmp.csv', 0, host_taxonomy_id, pc=self.pc, NoMSA=noMSA)
            data = Selenzy.updateScore(os.path.join(tmp_output_folder, '/tmp.csv'), score)
            val = json.loads(data.to_json())
            if 'Seq. ID' in val and len(val['Seq. ID'])>0:
                for ix in sorted(val['Seq. ID'], key=lambda z: int(z)):
                    uniprotID_score[val['Seq. ID'][ix]] = val['Score'][ix]
            else:
                raise ValueError
            return uniprotID_score


    def run(self,
            host_taxonomy_id=83333,
            pathway_id='rp_pathway',
            num_results=50,
            direction=0,
            noMSA=True,
            fp='RDK',
            rxntype='smarts',
            min_aa_length=100):
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
        for reac_id in self.getGroupsMembers(pathway_id):
            reac = self.model.getReaction(reac_id)
            brs_reac = self.readBRSYNTHAnnotation(reac.getAnnotation())
            if brs_reac['smiles']:
                try:
                    uniprotID_score = singleReactionRule(brs_reac['smiles'], host_taxonomy_id, num_results, direction, noMSA, fp, rxntype)
                    uniprotID_score_restricted = {}
                    for uniprot in uniprotID_score:
                        try:
                            if uniprot_aaLenght[uniprot]>int(min_aa_length):
                                uniprotID_score_restricted[uniprot] = uniprotID_score[uniprot]
                        except KeyError:
                            self.logger.warning('Cannot find the following UNIPROT '+str(uniprot)+' in uniprot_aaLenght')
                    xref = {'uniprot': [i for i in uniprotID_score_restricted]}
                    self.addUpdateMIRIAM(reac, 'reaction', xref)
                    self.addUpdateBRSynth(reac, 'selenzyme', uniprotID_score_restricted, None, False, True, True)
                except ValueError:
                    self.logger.warning('Problem with retreiving the selenzyme information for model '+str(self.model.getId()))
                    return False
            else:
                self.logger.warning('Cannot retreive the reaction rule of model '+str(self.model.getId()))
                return False
        return True
