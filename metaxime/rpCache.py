#!/usr/bin/env python3

import os
import csv
import logging
import numpy as np
import os
import pickle
import gzip
import urllib.request
import re
import tarfile
import shutil
import argparse
from ast import literal_eval
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import json

import objsize

__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = []
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Melchior du Lac"
__status__ = "Development"


#######################################################
################### rpCache  ##########################
#######################################################


class rpCache:
    """Class to generate the cache

    Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the the other steps. These should be called only when the files have changes

    """
    def __init__(self, rr_compounds_path=None, rr_rules_path=None, rr_rxn_recipes_path=None, fetch_input_files=False):
        """Cache constructor

        :param rr_compounds_path: Path to the compounds file
        :param rr_rules_path: Path to the rules file
        :param rr_rxn_recipes_path: Path to the reactions rules recipes
        :param fetch_input_files: Force download all the input cache file

        :type rr_compounds_path: str (Default: None)
        :type rr_rules_path: str (Default: None)
        :type rr_rxn_recipes_path: str (Default: None)
        :type fetch_input_files: bool (Default: False)

        :rtype: None
        :return: None
        """
        self.rr_compounds_path = rr_compounds_path
        self.rr_rules_path = rr_rules_path
        self.rr_rxn_recipes_path = rr_rxn_recipes_path
        #self.logger = logging.getLogger(__name__)
        self.logger = logging.getLogger(os.path.basename(__file__))
        self.logger.debug('Started instance of rpCache')
        #given by Thomas
        self.deprecatedCID_cid = {'MNXM162231': 'MNXM6',
                                  'MNXM84': 'MNXM15',
                                  'MNXM96410': 'MNXM14',
                                  'MNXM114062': 'MNXM3',
                                  'MNXM145523': 'MNXM57',
                                  'MNXM57425': 'MNXM9',
                                  'MNXM137': 'MNXM588022'}
        self.deprecatedRID_rid = {} #same structure as deprecatedCID_cid
        self.deprecatedComp_comp = {}
        #structure of the below looks like this: {'formula': 'C18H29O3', 'smiles': 'CCC=CCC1C(=O)CCC1CCCCCCCC(=O)O', 'inchi': 'InChI=1S/C18H30O3/c1-2-3-7-11-16-15(13-14-17(16)19)10-8-5-4-6-9-12-18(20)21/h3,7,15-16H,2,4-6,8-14H2,1H3,(H,20,21)', 'inchikey': 'BZXZFDKIRZBJEP-UHFFFAOYSA-N'}
        self.cid_strc = {} #cid_strc['MNXM1'] = {'formula': 'H', 'smiles': '[H+]', 'inchi': 'InChI=1S/p+1', 'inchikey': 'GPRLSGONYQIRFK-UHFFFAOYSA-N'}
        self.cid_xref = {} #cid_xref['MNXM287765'] = {'slm': ['000474284'], 'mnx': ['MNXM287765']}
        self.comp_xref = {}
        self.cid_name = {} #TODO: add valid names to the cid's
        #self.rid_name = {} #TODO: add valid reaction names
        self.chebi_cid = {} # for chebi_cid['88281']: 'MXM2323'
        self.inchikey_cid = {} # for key 'BZXZFDKIRZBJEP-UHFFFAOYSA-N': ['MNXM10', '10101']
        self.rr_reactions = {} # rr_reactions['RR-02-d2e7c5761b5a9b4b-04-F'] = {'MNXR139133': {'rule_id': 'RR-02-d2e7c5761b5a9b4b-04-F', 'rule_score': 0.3151075983206353, 'reac_id': 'MNXR139133', 'subs_id': 'MNXM89557', 'rel_direction': 1, 'left': {'MNXM89557': 1}, 'right': {'MNXM20': 1, 'MNXM722724': 1}}}
        self.rr_full_reactions = {} #rr_full_reactions['MNXR142257'] = {'left': {'MNXM4660': 1}, 'right': {'MNXM97172': 1}, 'direction': 0, 'main_left': ['MNXM4660'], 'main_right': ['MNXM97172']}
        self.dirname = os.path.dirname(os.path.abspath( __file__ ))
        if fetch_input_files:
            if not self._fetch_input_files():
                raise ValueError
        # cache
        if not os.path.isdir(os.path.join(self.dirname, 'cache')):
            os.mkdir(os.path.join(self.dirname, 'cache'))


    #####################################################
    ################# ERROR functions ###################
    #####################################################


    class Error(Exception):
        """Error function for the convertion of structures
        """
        pass


    class DepictionError(Error):
        """Error function for the convertion of structures
        """
        def __init__(self, message):
            """Constructor for the class

            :param message: The error handling message string

            :type message: str

            :rtype: None
            :return: None
            """
            #self.expression = expression
            self.message = message


    ##########################################################
    ################## Private Functions #####################
    ##########################################################


    ######## pubchem


    def _pubChemLimit(self):
        """Function to wait until the allowed number of requests can be made to pubChem

        No more than 5 requests per second.
        No more than 400 requests per minute.
        No longer than 300 second running time per minute.
        Requests exceeding limits are rejected (HTTP 503 error)

        :rtype: None
        :return: None
        """
        if self.pubchem_min_start==0.0:
            self.pubchem_min_start = time.time()
        self.pubchem_min_count += 1
        #### requests per minute ####
        if self.pubchem_min_count>=500 and time.time()-self.pubchem_min_start<=60.0:
            logger.warning('Reached 500 requests per minute for pubchem... waiting a minute')
            time.sleep(60.0)
            self.pubchem_min_start = time.time()
            self.pubchem_min_count = 0
        elif time.time()-self.pubchem_min_start>60.0:
            self.pubchem_min_start = time.time()
            self.pubchem_min_count = 0


    def _pubchemStrctSearch(self, strct, itype='inchi'):
        """Try to retreive the xref from an inchi structure using pubchem

        :param strct: The input structure
        :param itype: The type of input. Valid options: inchi, inchikey, smiles

        :type strct: str
        :type itype: str

        :rtype: dict
        :return: The resulting cross reference and structures
        """
        self._pubChemLimit()
        try:
            r = requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'+str(itype)+'/xrefs/SBURL/JSON', data={itype: strct})
            res_list = r.json()
        except json.decoder.JSONDecodeError:
            logger.warning('JSON decode error')
            return {}
        try:
            res_list = res_list['InformationList']['Information']
        except KeyError:
            logger.warning('pubchem JSON keyerror: '+str(res_list))
            return {}
        xref = {}
        if len(res_list)==1:
            _pubChemLimit()
            try:
                prop = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/property/IUPACName,InChI,InChIKey,CanonicalSMILES/JSON')
                prop_list = prop.json()
            except json.decoder.JSONDecodeError:
                logger.warning('JSON decode error')
                return {}
            try:
                name = prop_list['PropertyTable']['Properties'][0]['IUPACName']
                inchi = prop_list['PropertyTable']['Properties'][0]['InChI']
                inchikey = prop_list['PropertyTable']['Properties'][0]['InChIKey']
                smiles = prop_list['PropertyTable']['Properties'][0]['CanonicalSMILES']
            except KeyError:
                logger.warning('pubchem JSON keyerror: '+str(prop_list))
                return {}
            #TODO: need to determine how long cobra cannot handle this
            #TODO: determine if names that are too long is the problem and if not remove this part
            if len(name)>30:
                _pubChemLimit()
                try:
                    syn = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/synonyms/JSON')
                    syn_lst = syn.json()
                except json.decoder.JSONDecodeError:
                    logger.warning('pubchem JSON decode error')
                    return {}
                try:
                    syn_lst = syn_lst['InformationList']['Information'][0]['Synonym']
                    syn_lst = [x for x in syn_lst if not 'CHEBI' in x and not x.isupper()]
                    name = syn_lst[0] #need a better way instead of just the firs tone
                except KeyError:
                    logger.warning('pubchem JSON keyerror: '+str(syn.json()))
                    return {}
                except IndexError:
                    name = ''
            xref['pubchem'] = [str(res_list[0]['CID'])]
            for url in res_list[0]['SBURL']:
                if 'https://biocyc.org/compound?orgid=META&id=' in url:
                    if 'biocyc' not in xref:
                        xref['biocyc'] = []
                    xref['biocyc'].append(url.replace('https://biocyc.org/compound?orgid=META&id=', ''))
                if 'http://www.hmdb.ca/cidbolites/' in url:
                    if 'hmdb' not in xref:
                        xref['hmdb'] = []
                    xref['hmdb'].append(url.replace('http://www.hmdb.ca/cidbolites/', ''))
                if 'http://www.genome.jp/dbget-bin/www_bget?cpd:' in url:
                    if 'kegg_c' not in xref:
                        xref['kegg_c'] = []
                    xref['kegg_c'].append(url.replace('http://www.genome.jp/dbget-bin/www_bget?cpd:', ''))
                if 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' in url:
                    if 'chebi' not in xref:
                        xref['chebi'] = []
                    xref['chebi'].append(url.replace('http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:', ''))
        elif len(res_list)==0:
            logger.warning('Could not find results for: '+str(strct))
            return {}
        else:
            logger.warning('There are more than one result for '+str(strct)+'... Ignoring')
            return {}
        return {'name': name, 'inchi': inchi, 'inchikey': inchikey, 'smiles': smiles, 'xref': xref}


    def _fetch_input_files(self):
        """Private function to fetch the required data, parse them and generate the pickle

        :rtype: bool
        :return: Success or failure of the function
        """
        #################### make the local folders ############################
        # input_cache
        if not os.path.isdir(self.dirname+'/input_cache'):
            os.mkdir(self.dirname+'/input_cache')
        ####################### Fetch the input_cache files if necessary ######################
        # MNX 3.2
        url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'
        for in_file in ['reac_xref.tsv', 'chem_xref.tsv', 'chem_prop.tsv', 'comp_xref.tsv']:
            if not os.path.isfile(self.dirname+'/input_cache/'+in_file):
                urllib.request.urlretrieve(url+in_file, self.dirname+'/input_cache/'+in_file)
        # RetroRules
        if self.rr_compounds_path:
            shutil.copy(self.rr_compounds_path, os.path.join(self.dirname, 'input_cache', 'rr_compounds.tsv'))
        else:
            if not os.path.isfile(self.dirname+'/input_cache/rr_compounds.tsv'):
                urllib.request.urlretrieve('https://retrorules.org/dl/this/is/not/a/secret/path/rr02',
                                           self.dirname+'/input_cache/rr02_more_data.tar.gz')
                tar = tarfile.open(self.dirname+'/input_cache/rr02_more_data.tar.gz', 'r:gz')
                tar.extractall(self.dirname+'/input_cache/')
                tar.close()
                shutil.move(self.dirname+'/input_cache/rr02_more_data/compounds.tsv',
                            self.dirname+'/input_cache/rr_compounds.tsv')
                if not os.path.exists(self.dirname+'/input_cache/rxn_recipes.tsv'):
                    shutil.move(self.dirname+'/input_cache/rr02_more_data/rxn_recipes.tsv',
                                self.dirname+'/input_cache/')
                os.remove(self.dirname+'/input_cache/rr02_more_data.tar.gz')
                shutil.rmtree(self.dirname+'/input_cache/rr02_more_data')

        if self.rr_rxn_recipes_path:
            shutil.copy(self.rr_rxn_recipes_path, os.path.join(self.dirname, 'input_cache', 'rxn_recipes.tsv'))
        else:
            if not os.path.isfile(self.dirname+'/input_cache/rxn_recipes.tsv'):
                urllib.request.urlretrieve('https://retrorules.org/dl/this/is/not/a/secret/path/rr02',
                                           self.dirname+'/input_cache/rr02_more_data.tar.gz')
                tar = tarfile.open(self.dirname+'/input_cache/rr02_more_data.tar.gz', 'r:gz')
                tar.extractall(self.dirname+'/input_cache/')
                tar.close()
                if not os.path.exists(self.dirname+'/input_cache/rr_compounds.tsv'):
                    shutil.move(self.dirname+'/input_cache/rr02_more_data/compounds.tsv',
                                self.dirname+'/input_cache/rr_compounds.tsv')
                shutil.move(self.dirname+'/input_cache/rr02_more_data/rxn_recipes.tsv',
                            self.dirname+'/input_cache/')
                os.remove(self.dirname+'/input_cache/rr02_more_data.tar.gz')
                shutil.rmtree(self.dirname+'/input_cache/rr02_more_data')
        if self.rr_rules_path:
            shutil.copy(self.rr_rules_path, os.path.join(self.dirname, 'input_cache', 'rules_rall.tsv'))
        else:
            if not os.path.isfile(self.dirname+'/input_cache/rules_rall.tsv'):
                urllib.request.urlretrieve('https://retrorules.org/dl/preparsed/rr02/rp3/hs',
                                           self.dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz')
                tar = tarfile.open(self.dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz', 'r:gz')
                tar.extractall(self.dirname+'/input_cache/')
                tar.close()
                shutil.move(self.dirname+'/input_cache/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv', self.dirname+'/input_cache/rules_rall.tsv')
                os.remove(self.dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz')
                shutil.rmtree(self.dirname+'/input_cache/retrorules_rr02_rp3_hs')
        return True


    def _convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
        """Convert chemical depiction to others type of depictions

        Usage example:
         - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
         - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})

        :param idepic: Input string
        :param itype: The type of input
        :param otype: Type of output. Valid options: inchi, smiles, inchikey

        :type idepic: str
        :type itype: str
        :type otype: dict

        :rtype: dict
        :return: Dictionnary of results
        """
        # Import (if needed)
        self.logger.debug('input: '+str(idepic))
        self.logger.debug('itype: '+str(itype))
        self.logger.debug('otype: '+str(otype))
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise self.DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
        self.logger.debug('Sanitised the input')
        # Export
        odepic = dict()
        for item in otype:
            if item == 'smiles':
                odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
            elif item == 'inchi':
                odepic[item] = MolToInchi(rdmol)
            elif item == 'inchikey':
                odepic[item] = MolToInchiKey(rdmol)
            else:
                raise NotImplementedError('"{}" is not a valid output type'.format(otype))
        self.logger.debug('Exported the output')
        return odepic


    def _checkCIDdeprecated(self, cid):
        """Check if a given compound id is deprecated

        :param cid: Compound ID

        :type cid: str

        :rtype: str
        :return: The valid CID
        """
        if not self.deprecatedCID_cid:
            self.getDeprecatedCID()
        try:
            return self.deprecatedCID_cid[cid]
        except KeyError:
            return cid


    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkRIDdeprecated(self, rid):
        """Check if a given reaction id is deprecated

        :param rid: Reaction ID

        :type rid: str

        :rtype: str
        :return: The valid RID
        """
        if not self.deprecatedRID_rid:
            self.getDeprecatedRID()
        try:
            return self.deprecatedRID_rid[rid]
        except KeyError:
            return rid


    def _checkCompDeprecated(self, comp_id):
        """Check if the compartment id is deprecated

        :param comp_id: The compartment id

        :type comp_id: str

        :rtype: str
        :return: The compartment id
        """
        if not self.deprecatedComp_comp:
            self.getCompXref()
        try:
            return self.deprecatedComp_comp[comp_id]
        except KeyError:
            return comp_id


    #################################################################
    ################## Public functions #############################
    #################################################################


    ###################### Populate the cache #################################


    def populateCache(self):
        """Populate the cache by calling all the getters

        :rtype: None
        :return: None
        """
        self._fetch_input_files()
        a = self.getInchiKeyCID()
        a = self.getCIDxref()
        a = self.getCompXref()
        a = self.getCIDstrc()
        a = self.getDeprecatedCID()
        a = self.getDeprecatedRID()
        a = self.getFullReactions()
        a = self.getRRreactions()
        a = self.getChebiCID()
        a = self.getCIDname()


    def asDict(self, to_ret=None):
        """Return the whole set cache as a dictionnary
        """
        if not to_ret:
            to_ret = {}
        to_ret['cid_strc'] = self.cid_strc
        to_ret['deprecatedCID_cid'] = self.deprecatedCID_cid
        to_ret['deprecatedRID_rid'] = self.deprecatedRID_rid
        to_ret['cid_xref'] = self.cid_xref
        to_ret['comp_xref'] = self.comp_xref
        to_ret['rr_full_reactions'] = self.rr_full_reactions
        to_ret['chebi_cid'] = self.chebi_cid
        to_ret['inchikey_cid'] = self.inchikey_cid
        to_ret['cid_name'] = self.cid_name
        return to_ret


    def setFromDict(self, cache_dict):
        """Set the current cache from a dict input
        """
        try:
            self.cid_strc = cache_dict['cid_strc']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for cid_strc')
        try:
            self.deprecatedCID_cid = cache_dict['deprecatedCID_cid']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for deprecatedCID_cid')
        try:
            self.deprecatedRID_rid = cache_dict['deprecatedRID_rid']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for deprecatedRID_rid')
        try:
            self.cid_xref = cache_dict['cid_xref']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for cid_xref')
        try:
            self.comp_xref = cache_dict['comp_xref']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for comp_xref')
        try:
            self.rr_full_reactions = cache_dict['rr_full_reactions']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for rr_full_reactions')
        try:
            self.chebi_cid = cache_dict['chebi_cid']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for chebi_cid')
        try:
            self.inchikey_cid = cache_dict['inchikey_cid']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for inchikey_cid')
        try:
            self.cid_name = cache_dict['cid_name']
        except KeyError:
            self.logger.warning('There does not seem to be an entry for cid_name')
        return True


    def getSizeCache(self):
        size_cache = 0
        size_cache += objsize.get_deep_size(self.cid_strc)
        size_cache += objsize.get_deep_size(self.deprecatedCID_cid)
        size_cache += objsize.get_deep_size(self.deprecatedRID_rid)
        size_cache += objsize.get_deep_size(self.cid_xref)
        size_cache += objsize.get_deep_size(self.comp_xref)
        size_cache += objsize.get_deep_size(self.rr_full_reactions)
        size_cache += objsize.get_deep_size(self.chebi_cid)
        size_cache += objsize.get_deep_size(self.inchikey_cid)
        size_cache += objsize.get_deep_size(self.cid_name)
        return size_cache


    ##################### Individual getters for the cache ####################


    ### MNX
    def getCIDstrc(self):
        """Getter of the dictionary describing the structure of compound id's

        :rtype: dict
        :return: The dictionnary structure of compound
        """
        #this one hos both mnx and RR
        picklename = 'cid_strc.pickle.gz'
        filename = 'chem_prop.tsv'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            self.MNXstrc(os.path.join(self.dirname, 'input_cache', filename))
            self.retroRulesStrc(self.dirname+'/input_cache/rr_compounds.tsv')
            pickle.dump(self.cid_strc, gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.cid_strc = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.cid_strc
        #return pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))


    def getDeprecatedCID(self):
        """Getter of the dictionary describing the deprecated compound id to compound id

        :rtype: dict
        :return: The dictionnary structure of deprecated compounds
        """
        picklename = 'deprecatedCID_cid.pickle.gz'
        filename = 'chem_xref.tsv'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            self.deprecatedMNXM(os.path.join(self.dirname, 'input_cache', filename))
            pickle.dump(self.deprecatedCID_cid, gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.deprecatedCID_cid = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.deprecatedCID_cid
        #return pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))


    def getDeprecatedRID(self):
        """Getter of the dictionary describing the deprecated reaction id to reaction id

        :rtype: dict
        :return: The dictionnary structure of deprecated reactions
        """
        picklename = 'deprecatedRID_rid.pickle.gz'
        filename = 'reac_xref.tsv'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            self.deprecatedMNXR(os.path.join(self.dirname, 'input_cache', filename))
            pickle.dump(self.deprecatedRID_rid, gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.deprecatedRID_rid = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.deprecatedRID_rid
        #return pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))


    def getCIDxref(self):
        """Getter of the dictionary describing the cross-reference of compounds

        :rtype: dict
        :return: The dictionnary of cross-reference
        """
        picklename = 'cid_xref.pickle.gz'
        filename = 'chem_xref.tsv'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            self.mnxXref(os.path.join(self.dirname, 'input_cache', filename))
            pickle.dump(self.cid_xref,
                        gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.cid_xref = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.cid_xref
        #return pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))


    def getCompXref(self):
        """Getter of the dictionary describing the cross-reference of the compartment

        :rtype: dict
        :return: The dictionnary of cross-reference
        """
        picklename1 = 'xref_comp.pickle.gz'
        picklename2 = 'comp_xref.pickle.gz'
        filename = 'comp_xref.tsv'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename1)) or not os.path.isfile(os.path.join(self.dirname, 'cache', picklename2)):
            self.mnxCompXref(self.dirname+'/input_cache/comp_xref.tsv')
            pickle.dump(self.deprecatedComp_comp, gzip.open(os.path.join(self.dirname, 'cache', picklename1),'wb'))
            pickle.dump(self.comp_xref, gzip.open(os.path.join(self.dirname, 'cache', picklename2), 'wb'))
        self.deprecatedComp_comp = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename1), 'rb'))
        self.comp_xref = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename2), 'rb'))
        return self.deprecatedComp_comp, self.comp_xref


    def getFullReactions(self):
        """Getter of the dictionary describing the full reactions from reaction rules

        :rtype: dict
        :return: The dictionnary of full reactions
        """
        picklename = 'rr_full_reactions.pickle.gz'
        filename = 'rxn_recipes.tsv'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            self.retroRulesFullReac(os.path.join(self.dirname, 'input_cache', filename))
            pickle.dump(self.rr_full_reactions,
                        gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.rr_full_reactions = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.rr_full_reactions


    def getRRreactions(self):
        """Getter of the dictionary describing the reaction rules

        :rtype: dict
        :return: The dictionnary of reaction rules
        """
        picklename = 'rr_reactions.pickle.gz'
        filename = 'rules_rall.tsv'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            self.retroReactions(os.path.join(self.dirname, 'input_cache', filename))
            pickle.dump(self.rr_reactions,
                        gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.rr_reactions = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.rr_reactions


    ### intermidiate generate, WARNING: needs to be added at the end after the other files
    def getChebiCID(self):
        """Getter of the dictionnary describing the Chebi id to the compound id

        :rtype: dict
        :return: The dictionnary of Chebi CID
        """
        picklename = 'chebi_cid.pickle.gz'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            if not self.cid_xref:
                self.getCIDxref()
            for cid in self.cid_xref:
                if 'chebi' in self.cid_xref[cid]:
                    for c in self.cid_xref[cid]['chebi']:
                        self.chebi_cid[c] = cid
            pickle.dump(self.chebi_cid,
                        gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.chebi_cid = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.chebi_cid


    def getInchiKeyCID(self):
        """Getter of the dictionnary describing the inchikey to the compound id

        :rtype: dict
        :return: The dictionnary of inchikey to cid
        """
        picklename = 'inchikey_cid.pickle.gz'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            if not self.cid_strc:
                self.getCIDstrc()
            # open the already calculated (normally) mnxm_strc.pickle.gz
            for cid in self.cid_strc:
                if not self.cid_strc[cid]['inchikey'] in self.inchikey_cid:
                    self.inchikey_cid[self.cid_strc[cid]['inchikey']] = []
                if not cid in self.inchikey_cid[self.cid_strc[cid]['inchikey']]:
                    self.inchikey_cid[self.cid_strc[cid]['inchikey']].append(cid)
            pickle.dump(self.inchikey_cid, gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.inchikey_cid = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.inchikey_cid


    def getCIDname(self):
        """Getter of the dictionnary describing the name of the compounds

        :rtype: dict
        :return: The dictionnary of compound name
        """
        picklename = 'cid_name.pickle.gz'
        if not os.path.isfile(os.path.join(self.dirname, 'cache', picklename)):
            pickle.dump(self.cid_name,
                        gzip.open(os.path.join(self.dirname, 'cache', picklename), 'wb'))
        self.cid_name = pickle.load(gzip.open(os.path.join(self.dirname, 'cache', picklename), 'rb'))
        return self.cid_name


    #TODO
    #def getRIDxref(self):


    ################## RetroRules parsers ####################################


    def retroRulesStrc(self, rr_compounds_path):
        """Parse the Reaction Rules

        :param rr_compounds_path: The reaction rules compounds

        :type cc_compounds_path: str

        :rtype: None
        :return: None
        """
        for row in csv.DictReader(open(rr_compounds_path), delimiter='\t'):
            tmp = {'formula':  None,
                   'smiles': None,
                   'inchi': row['inchi'],
                   'inchikey': row['inchikey']}
            try:
                resConv = self._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except self.DepictionError as e:
                self.logger.warning('Could not convert some of the structures: '+str(tmp))
                self.logger.warning(e)
            if not self._checkCIDdeprecated(row['cid']) in self.cid_strc:
                self.cid_strc[self._checkCIDdeprecated(row['cid'])] = tmp
            else:
                if resConv['smiles'] and not self.cid_strc[self._checkCIDdeprecated(row['cid'])]['smiles']:
                    self.cid_strc[self._checkCIDdeprecated(row['cid'])]['smiles'] = resConv['smiles']


    #NOTE: we take care of the fact that a given reaction rule may have multiple reactions associated with them
    def retroReactions(self, rules_rall_path):
        """Function to parse the rules_rall.tsv from RetroRules

        Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction

        :param rules_rall_path: Path to the reaction rules

        :type rules_rall_path: str

        :rtype: bool
        :return: Success or failure of the function
        """
        try:
            for row in csv.DictReader(open(rules_rall_path), delimiter='\t'):
                if not type(row)==dict:
                    row = dict(row)
                #NOTE: as of now all the rules are generated using MNX
                #but it may be that other db are used, we are handling this case
                #WARNING: can have multiple products so need to seperate them
                products = {}
                for i in row['Product_IDs'].split('.'):
                    cid = self._checkCIDdeprecated(i)
                    if not cid in products:
                        products[cid] = 1
                    else:
                        products[cid] += 1
                try:
                    #WARNING: one reaction rule can have multiple reactions associated with them
                    #To change when you can set subpaths from the mutliple numbers of
                    #we assume that the reaction rule has multiple unique reactions associated
                    if row['# Rule_ID'] not in self.rr_reactions:
                        self.rr_reactions[row['# Rule_ID']] = {}
                    if row['# Rule_ID'] in self.rr_reactions[row['# Rule_ID']]:
                        self.logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                    self.rr_reactions[row['# Rule_ID']][self._checkRIDdeprecated(row['Reaction_ID'])] = {'rule_id': row['# Rule_ID'],
                                                                                                         'rule_score': float(row['Score_normalized']),
                                                                                                         'reac_id': self._checkRIDdeprecated(row['Reaction_ID']),
                                                                                                         'subs_id': self._checkCIDdeprecated(row['Substrate_ID']),
                                                                                                         'rel_direction': int(row['Rule_relative_direction']),
                                                                                                         'left': {self._checkCIDdeprecated(row['Substrate_ID']): 1},
                                                                                                         'right': products}
                except ValueError:
                    self.logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    self.logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
        except FileNotFoundError as e:
                self.logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
        return True


    def retroRulesFullReac(self, rxn_recipes_path):
        """Generate complete reactions from the rxn_recipes.tsv from RetroRules

        These are the compplete reactions from which the reaction rules are generated from. This is used to reconstruct the full reactions from monocomponent reactions

        :param rxn_recipes_path: Path to the reaction recipies

        :type rxn_recipes_path: str

        :rtype: None
        :return: None
        """
        #### for character matching that are returned
        DEFAULT_STOICHIO_RESCUE = {'4n': 4, '3n': 3, "2n": 2, 'n': 1,
                                   '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
                                   'N': 1, 'm': 1, 'q': 1,
                                   '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
                                   '0.02': 1, '0.2': 1,
                                   '(n-1)': 0, '(n-2)': -1}
        try:
            for row in csv.DictReader(open(rxn_recipes_path), delimiter='\t'):
                tmp = {} # makes sure that if theres an error its not added
                #parse the reaction equation
                if not len(row['Equation'].split('='))==2:
                    self.logger.warning('There should never be more or less than a left and right of an euation')
                    self.logger.warnin(row['Equation'])
                    continue
                ######### LEFT ######
                #### MNX id
                tmp['left'] = {}
                for spe in re.findall('(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[0]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['left'][self._checkCIDdeprecated(spe[1])] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['left'][self._checkCIDdeprecated(spe[1])] = int(spe[0])
                        except ValueError:
                            self.logger.warning('Cannot convert '+str(spe[0]))
                            continue
                ####### RIGHT #####
                ####  MNX id
                tmp['right'] = {}
                for spe in re.findall('(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[1]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['right'][self._checkCIDdeprecated(spe[1])] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['right'][self._checkCIDdeprecated(spe[1])] = int(spe[0])
                        except ValueError:
                            self.logger.warning('Cannot convert '+str(spe[0]))
                            continue
                ####### DIRECTION ######
                try:
                    tmp['direction'] = int(row['Direction'])
                except ValueError:
                    self.logger.error('Cannot convert '+str(row['Direction'])+' to int')
                    continue
                ### add the others
                tmp['main_left'] = row['Main_left'].split(',')
                tmp['main_right'] = row['Main_right'].split(',')
                self.rr_full_reactions[self._checkRIDdeprecated(row['#Reaction_ID'])] = tmp
        except FileNotFoundError:
            self.logger.error('Cannot find file: '+str(path))


    ################## MNX parsers ###################################


    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def mnxCompXref(self, comp_xref_path):
        """Function to parse the comp_xref.tsv file of MetanetX

        Generate a dictionnary of compartments id's (MNX) to other database id's

        :param comp_xref_path: Path to the compartment xref

        :type comp_xref_path: str

        :rtype: None
        :return: None
        """
        try:
            with open(comp_xref_path) as f:
                c = csv.reader(f, delimiter='\t')
                #not_recognised = []
                for row in c:
                    #cid = row[0].split(':')
                    if not row[0][0]=='#':
                        #collect the info
                        mnxc = row[1]
                        if len(row[0].split(':'))==1:
                            dbName = 'mnx'
                            dbCompId = row[0]
                        else:
                            dbName = row[0].split(':')[0]
                            dbCompId = ''.join(row[0].split(':')[1:])
                            dbCompId = dbCompId.lower()
                        if dbName=='deprecated':
                            dbName = 'mnx'
                        #create the dicts
                        if not mnxc in self.comp_xref:
                            self.comp_xref[mnxc] = {}
                        if not dbName in self.comp_xref[mnxc]:
                            self.comp_xref[mnxc][dbName] = []
                        if not dbCompId in self.comp_xref[mnxc][dbName]:
                            self.comp_xref[mnxc][dbName].append(dbCompId)
                        #create the reverse dict
                        if not dbCompId in self.deprecatedComp_comp:
                            self.deprecatedComp_comp[dbCompId] = mnxc
        except FileNotFoundError:
            self.logger.error('comp_xref file not found')
            return {}



    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
    def deprecatedMNXM(self, chem_xref_path):
        """Function to parse the chem_xref.tsv file of MetanetX

        Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's. This can include more than one old id per new one and thus returns a dictionnary. Private function

        :param chem_xref_path: The path to the MetaNetX chemical list

        :type chem_xref_path: str

        :rtype: None
        :return: None
        """
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        self.deprecatedCID_cid[mnx[1]] = row[1]
            try:
                self.deprecatedCID_cid['MNXM01'] = 'MNXM1'
            except KeyError:
                pass


    def deprecatedMNXR(self, reac_xref_path):
        """Function to parse the reac_xref.tsv file of MetanetX

        Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's. This can include more than one old id per new one and thus returns a dictionnary. Private function

        :param reac_xref_path: The path to the MetaNetX reaction list

        :type reac_xref_path: str

        :rtype: None
        :return: None
        """
        with open(reac_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        self.deprecatedRID_rid[mnx[1]] = row[1]


    def MNXstrc(self, chem_prop_path):
        """Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules

        Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components

        :param chem_prop_path: The path to the MetaNetX chemical structure list

        :type chem_prop_path: str

        :rtype: None
        :return: None
        """
        with open(chem_prop_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    self.logger.debug('--- Parsing '+str(row[0])+' ----') 
                    mnxm = self._checkCIDdeprecated(row[0])
                    tmp = {'formula':  row[2],
                           'name': row[1],
                           'smiles': row[6],
                           'inchi': row[5],
                           'inchikey': row[8]}
                    for i in tmp:
                        if tmp[i]=='' or tmp[i]=='NA':
                            tmp[i] = None
                    if not mnxm in self.cid_name and tmp['name']:
                        self.cid_name[mnxm] = tmp['name']
                    if mnxm in self.cid_strc:
                        self.cid_strc[mnxm]['formula'] = row[2]
                        if not self.cid_strc[mnxm]['smiles'] and tmp['smiles']:
                            self.cid_strc[mnxm]['smiles'] = tmp['smiles']
                        if not self.cid_strc[mnxm]['inchikey'] and tmp['inchikey']:
                            self.cid_strc[mnxm]['inchikey'] = tmp['inchikey']
                    else:
                        #check to see if the inchikey is valid or not
                        otype = set({})
                        if not tmp['inchikey']:
                            otype.add('inchikey')
                        if not tmp['smiles']:
                            otype.add('smiles')
                        if not tmp['inchi']:
                            otype.add('inchi')
                        itype = ''
                        if tmp['inchi']:
                            itype = 'inchi'
                        elif tmp['smiles']:
                            itype = 'smiles'
                        else:
                            self.logger.warning('No valid entry for the convert_depiction function')
                            continue
                        if otype:
                            try:
                                self.logger.debug('Converting using '+str(tmp[itype]))
                                resConv = self._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                                for i in resConv:
                                    tmp[i] = resConv[i]
                            except self.DepictionError as e:
                                self.logger.warning('Could not convert some of the structures: '+str(tmp))
                                self.logger.warning(e)
                        self.cid_strc[mnxm] = tmp


    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
    def mnxXref(self, chem_xref_path):
        """Function to parse the chem_xref.tsv file of MetanetX

        Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's

        :param chem_xref_path: The path to the MetaNetX chemical cross reference list

        :type chem_xref_path: str

        :rtype: None
        :return: None
        """
        with open(chem_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = self._checkCIDdeprecated(row[1])
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    #mnx
                    if not mnx in self.cid_xref:
                        self.cid_xref[mnx] = {}
                    if not dbName in self.cid_xref[mnx]:
                        self.cid_xref[mnx][dbName] = []
                    if not dbId in self.cid_xref[mnx][dbName]:
                        self.cid_xref[mnx][dbName].append(dbId)
                    ### DB ###
                    if not dbName in self.cid_xref:
                        self.cid_xref[dbName] = {}
                    if not dbId in self.cid_xref[dbName]:
                        self.cid_xref[dbName][dbId] = mnx


    ##################### QUERY #####################

    def queryRRreactions(self, rr):
        """Query the structure information of a chemical species

        :param cid: A chemical id

        :type cid: str

        :rtype: dict
        :return: Dictionnary of the known chemical structures
        """
        if not self.rr_reactions:
            self.getRRreactions()
        try:
            return self.rr_reactions[rr]
        except KeyError:
            self.logger.warning('Cache does have any information for the following reaction rule id: '+str(rr))
            return {}

    def queryCIDstrc(self, cid):
        """Query the structure information of a chemical species

        :param cid: A chemical id

        :type cid: str

        :rtype: dict
        :return: Dictionnary of the known chemical structures
        """
        if not self.cid_strc:
            self.getCIDstrc()
        if not self.deprecatedCID_cid:
            self.getDeprecatedCID()
        try:
            return self.cid_strc[self._checkCIDdeprecated(cid)]
        except KeyError:
            self.logger.warning('Cache does have structure information for: '+str(cid))
            return {}


    def queryChebiCID(self, chebi):
        """Query the chemical id from a chebi id

        :param chebi: A chebi chemical id

        :type chebi: str

        :rtype: str
        :return: The chemical id
        """
        if not self.chebi_cid:
            self.getChebiCID()
        if not self.deprecatedCID_cid:
            self.getDeprecatedCID()
        try:
            return self._checkCIDdeprecated(self.chebi_cid[chebi])
        except KeyError:
            self.logger.warning('Cache does not have CHEBI entry: '+str(chebi))
            return None


    def queryCompXref(self, comp_id):
        """Query the compartment id

        :param comp_id: A compartment id

        :type comp_id: str

        :rtype: dict
        :return: The list of compartment id's
        """
        if not self.comp_xref:
            self.getCompXref()
        try:
            return self.comp_xref[comp_id]
        except KeyError:
            self.logger.warning('Cache does not have the compartment entry: '+str(comp_id))
            return {}


    def queryCIDxref(self, cid):
        """Query the chemical species cross-reference

        :param cid: A chemical species id

        :type cid: str

        :rtype: dict
        :return: The cross-reference of that species
        """
        if not self.cid_xref:
            self.getCIDxref()
        if not self.deprecatedCID_cid:
            self.getDeprecatedCID()
        try:
            return self.cid_xref[self._checkCIDdeprecated(cid)]
        except KeyError:
            self.logger.warning('Cache does not have the chemical species: '+str(cid))
            return {}


    def queryCIDname(self, cid):
        """Query the name of a chemical species

        :param cid: A chemical id

        :type cid: str

        :rtype: str
        :return: Common given name
        """
        if not self.cid_name:
            self.getCIDname()
        if not self.deprecatedCID_cid:
            self.getDeprecatedCID()
        try:
            return self.cid_name[self._checkCIDdeprecated(cid)]
        except KeyError:
            self.logger.warning('Cache does have structure information for: '+str(cid))
            return {}


    def queryInchiKeyCID(self, inchikey):
        """Query the chemical id for the inchikey

        :param inchikey: An inchikey

        :type inchikey: str

        :rtype: str
        :return: Chemical ID
        """
        if not self.inchikey_cid:
            self.getInchiKeyCID()
        try:
            return self._checkCIDdeprecated([i for i in self.inchikey_cid[inchikey] if i[:3]=='MNX'][0])
        except KeyError:
            self.logger.warning('Cache does have the following CID for the inchikey: '+str(inchikey))
            return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Pass different Rules files')
    parser.add_argument('-rr_compounds_path', type=str, default=None)
    parser.add_argument('-rr_rules_path', type=str, default=None)
    parser.add_argument('-rr_rxn_recipes_path', type=str, default=None)
    parser.add_argument('-fetch_input_files', type=str, default='True')
    params = parser.parse_args()
    if params.fetch_input_files==True or params.fetch_input_files=='True' or params.fetch_input_files=='true':
        fetch_input_files = True
    elif params.fetch_input_files==False or params.fetch_input_files=='False' or params.fetch_input_files=='false':
        fetch_input_files = False
    else:
        logging.error('Cannot interpret '+str(params.fetch_input_files))
        exit(1)
    if params.rr_compounds_path:
        if not os.path.exists(params.rr_compounds_path):
            logging.error('The file: '+str(params.rr_compounds_path)+' does not exist')
            exit(1)
    if params.rr_rules_path:
        if not os.path.exists(params.rr_rules_path):
            logging.error('The file: '+str(params.rr_rules_path)+' does not exist')
            exit(1)
    if params.rr_rxn_recipes_path:
        if not os.path.exists(params.rr_rxn_recipes_path):
            logging.error('The file: '+str(params.rr_rxn_recipes_path)+' does not exist')
            exit(1)
    rpcache = rpCache(params.rr_compounds_path, params.rr_rules_path, params.rr_rxn_recipes_path, fetch_input_files)
    rpcache.populateCache()
