import os
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
import csv
import logging
import os
import pickle
import gzip
import urllib.request
import re
import tarfile
import shutil
import argparse

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

'''
logging.disable(logging.INFO)
logging.disable(logging.WARNING)
'''


#######################################################
################### rpCache  ##########################
#######################################################


## Class to generate the cache
#
# Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the
#the other steps. These should be called only when the files have changes
class rpCache:
    ## Cache constructor
    #
    # @param self The object pointer
    # @param inputPath The path to the folder that contains all the input/output files required
    def __init__(self, rr_compounds_path=None, rr_rules_path=None, rr_rxn_recipes_path=None):
        self.rr_compounds_path = rr_compounds_path
        self.rr_rules_path = rr_rules_path
        self.rr_rxn_recipes_path = rr_rxn_recipes_path
        #given by Thomas
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpCache')
        self.deprecatedCID_cid = {'MNXM162231': 'MNXM6',
                                  'MNXM84': 'MNXM15',
                                  'MNXM96410': 'MNXM14',
                                  'MNXM114062': 'MNXM3',
                                  'MNXM145523': 'MNXM57',
                                  'MNXM57425': 'MNXM9',
                                  'MNXM137': 'MNXM588022'}
        self.deprecatedRID_rid = {} #same structure as deprecatedCID_cid
        #structure of the below looks like this: {'formula': 'C18H29O3', 'smiles': 'CCC=CCC1C(=O)CCC1CCCCCCCC(=O)O', 'inchi': 'InChI=1S/C18H30O3/c1-2-3-7-11-16-15(13-14-17(16)19)10-8-5-4-6-9-12-18(20)21/h3,7,15-16H,2,4-6,8-14H2,1H3,(H,20,21)', 'inchikey': 'BZXZFDKIRZBJEP-UHFFFAOYSA-N'}
        self.cid_strc = {} #cid_strc['MNXM1'] = {'formula': 'H', 'smiles': '[H+]', 'inchi': 'InChI=1S/p+1', 'inchikey': 'GPRLSGONYQIRFK-UHFFFAOYSA-N'}
        self.cid_xref = {} #cid_xref['MNXM287765'] = {'slm': ['000474284'], 'mnx': ['MNXM287765']}
        self.comp_xref = {}
        self.xref_comp = {}
        #self.cid_name = {} #TODO: add valid names to the cid's
        #self.rid_name = {} #TODO: add valid reaction names
        self.chebi_cid = {} # for chebi_cid['88281']: 'MXM2323'
        self.inchikey_cid = {} # for key 'BZXZFDKIRZBJEP-UHFFFAOYSA-N': ['MNXM10', '10101']
        self.rr_reactions = {} # rr_reactions['RR-02-d2e7c5761b5a9b4b-04-F'] = {'MNXR139133': {'rule_id': 'RR-02-d2e7c5761b5a9b4b-04-F', 'rule_score': 0.3151075983206353, 'reac_id': 'MNXR139133', 'subs_id': 'MNXM89557', 'rel_direction': 1, 'left': {'MNXM89557': 1}, 'right': {'MNXM20': 1, 'MNXM722724': 1}}}
        self.rr_full_reactions = {} #rr_full_reactions['MNXR142257'] = {'left': {'MNXM4660': 1}, 'right': {'MNXM97172': 1}, 'direction': 0, 'main_left': ['MNXM4660'], 'main_right': ['MNXM97172']}
        if not self._loadCache():
            raise ValueError

    #####################################################
    ################# ERROR functions ###################
    #####################################################

    ## Error function for the convertion of structures
    #
    class Error(Exception):
        pass


    ## Error function for the convertion of structures
    #
    class DepictionError(Error):
        def __init__(self, message):
            #self.expression = expression
            self.message = message


    ##########################################################
    ################## Private Functions #####################
    ##########################################################

    ## Private function to fetch the required data, parse them and generate the pickle
    #
    #  Opens the previously generated cache to the object memory
    #
    # @param The oject pointer
    # @return Boolean detemining the success of the function or not
    def _loadCache(self, fetchInputFiles=False):
        dirname = os.path.dirname(os.path.abspath( __file__ ))
        #################### make the local folders ############################
        # input_cache
        if not os.path.isdir(dirname+'/input_cache'):
            os.mkdir(dirname+'/input_cache')
        # cache
        if not os.path.isdir(dirname+'/cache'):
            os.mkdir(dirname+'/cache')

        ####################### Fetch the input_cache files if necessary ######################

        # MNX
        url = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/'
        for in_file in ['reac_xref.tsv', 'chem_xref.tsv', 'chem_prop.tsv', 'comp_xref.tsv']:
            if not os.path.isfile(dirname+'/input_cache/'+in_file) or fetchInputFiles:
                urllib.request.urlretrieve(url+in_file, dirname+'/input_cache/'+in_file)

        # RetroRules
        if self.rr_compounds_path:
            shutil.copy(self.rr_compounds_path, os.path.join(dirname, 'input_cache', 'rr_compounds.tsv'))
        else:
            if not os.path.isfile(dirname+'/input_cache/rr_compounds.tsv') or fetchInputFiles:
                urllib.request.urlretrieve('https://retrorules.org/dl/this/is/not/a/secret/path/rr02',
                                           dirname+'/input_cache/rr02_more_data.tar.gz')
                tar = tarfile.open(dirname+'/input_cache/rr02_more_data.tar.gz', 'r:gz')
                tar.extractall(dirname+'/input_cache/')
                tar.close()
                shutil.move(dirname+'/input_cache/rr02_more_data/compounds.tsv',
                            dirname+'/input_cache/rr_compounds.tsv')
                if not os.path.exists(dirname+'/input_cache/rxn_recipes.tsv'):
                    shutil.move(dirname+'/input_cache/rr02_more_data/rxn_recipes.tsv',
                                dirname+'/input_cache/')
                os.remove(dirname+'/input_cache/rr02_more_data.tar.gz')
                shutil.rmtree(dirname+'/input_cache/rr02_more_data')

        if self.rr_rxn_recipes_path:
            shutil.copy(self.rr_rxn_recipes_path, os.path.join(dirname, 'input_cache', 'rxn_recipes.tsv'))
        else:
            if not os.path.isfile(dirname+'/input_cache/rxn_recipes.tsv') or fetchInputFiles:
                urllib.request.urlretrieve('https://retrorules.org/dl/this/is/not/a/secret/path/rr02',
                                           dirname+'/input_cache/rr02_more_data.tar.gz')
                tar = tarfile.open(dirname+'/input_cache/rr02_more_data.tar.gz', 'r:gz')
                tar.extractall(dirname+'/input_cache/')
                tar.close()
                if not os.path.exists(dirname+'/input_cache/rr_compounds.tsv'):
                    shutil.move(dirname+'/input_cache/rr02_more_data/compounds.tsv',
                                dirname+'/input_cache/rr_compounds.tsv')
                shutil.move(dirname+'/input_cache/rr02_more_data/rxn_recipes.tsv',
                            dirname+'/input_cache/')
                os.remove(dirname+'/input_cache/rr02_more_data.tar.gz')
                shutil.rmtree(dirname+'/input_cache/rr02_more_data')

        if self.rr_rules_path:
            shutil.copy(self.rr_rules_path, os.path.join(dirname, 'input_cache', 'rules_rall.tsv'))
        else:
            if not os.path.isfile(dirname+'/input_cache/rules_rall.tsv') or fetchInputFiles:
                urllib.request.urlretrieve('https://retrorules.org/dl/preparsed/rr02/rp3/hs',
                                           dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz')
                tar = tarfile.open(dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz', 'r:gz')
                tar.extractall(dirname+'/input_cache/')
                tar.close()
                shutil.move(dirname+'/input_cache/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv', dirname+'/input_cache/rules_rall.tsv')
                os.remove(dirname+'/input_cache/retrorules_rr02_rp3_hs.tar.gz')
                shutil.rmtree(dirname+'/input_cache/retrorules_rr02_rp3_hs')


        ###################### Populate the cache #################################

        ###### RetroRules #####
        picklename = 'rr_full_reactions.pickle.gz'
        filename = 'rxn_recipes.tsv'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename)):
            self.retroRulesFullReac(os.path.join(dirname, 'input_cache', filename))
            pickle.dump(self.rr_full_reactions,
                        gzip.open(os.path.join(dirname, 'cache', picklename), 'wb'))
        self.rr_full_reactions = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename), 'rb'))

        picklename = 'rr_reactions.pickle.gz'
        filename = 'rules_rall.tsv'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename)):
            self.retroReactions(os.path.join(dirname, 'input_cache', filename))
            pickle.dump(self.rr_reactions,
                        gzip.open(os.path.join(dirname, 'cache', picklename), 'wb'))
        self.rr_reactions = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename), 'rb'))

        ### MNX

        #this one hos both mnx and RR
        picklename = 'cid_strc.pickle.gz'
        filename = 'chem_prop.tsv'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename)):
            self.MNXstrc(os.path.join(dirname, 'input_cache', filename))
            self.retroRulesStrc(dirname+'/input_cache/rr_compounds.tsv')
            pickle.dump(self.cid_strc, gzip.open(os.path.join(dirname, 'cache', picklename), 'wb'))
        self.cid_strc = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename), 'rb'))

        picklename = 'deprecatedCID_cid.pickle.gz'
        filename = 'chem_xref.tsv'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename)):
            self.deprecatedMNXM(os.path.join(dirname, 'input_cache', filename))
            pickle.dump(self.deprecatedCID_cid, gzip.open(os.path.join(dirname, 'cache', picklename), 'wb'))
        self.deprecatedCID_cid = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename), 'rb'))

        picklename = 'deprecatedRID_rid.pickle.gz'
        filename = 'reac_xref.tsv'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename)):
            self.deprecatedMNXR(os.path.join(dirname, 'input_cache', filename))
            pickle.dump(self.deprecatedRID_rid, gzip.open(os.path.join(dirname, 'cache', picklename), 'wb'))
        self.deprecatedRID_rid = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename), 'rb'))

        picklename = 'cid_xref.pickle.gz'
        filename = 'chem_xref.tsv'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename)):
            self.mnxXref(os.path.join(dirname, 'input_cache', filename))
            pickle.dump(self.cid_xref,
                        gzip.open(os.path.join(dirname, 'cache', picklename), 'wb'))
        self.cid_xref = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename), 'rb'))

        picklename1 = 'xref_comp.pickle.gz'
        picklename2 = 'comp_xref.pickle.gz'
        filename = 'comp_xref.tsv'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename1)) or not os.path.isfile(os.path.join(dirname, 'cache', picklename2)):
            self.mnxCompXref(dirname+'/input_cache/comp_xref.tsv')
            pickle.dump(self.xref_comp, gzip.open(os.path.join(dirname, 'cache', picklename1),'wb'))
            pickle.dump(self.comp_xref, gzip.open(os.path.join(dirname, 'cache', picklename2), 'wb'))
        self.xref_comp = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename1), 'rb'))
        self.comp_xref = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename2), 'rb'))

        ### intermidiate generate

        picklename = 'chebi_cid.pickle.gz'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename)):
            self._chebiXref()
            pickle.dump(self.chebi_cid,
                        gzip.open(os.path.join(dirname, 'cache', picklename), 'wb'))
        self.chebi_cid = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename), 'rb'))

        picklename = 'inchikey_cid.pickle.gz'
        if not os.path.isfile(os.path.join(dirname, 'cache', picklename)):
            # open the already calculated (normally) mnxm_strc.pickle.gz
            self._inchikeyCID()
            pickle.dump(self.inchikey_cid, gzip.open(os.path.join(dirname, 'cache', picklename), 'wb'))
        self.inchikey_cid = pickle.load(gzip.open(os.path.join(dirname, 'cache', picklename), 'rb'))
        return True


    #####################################################
    ############## Private Function #####################
    #####################################################


    ## Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
    #  @param self The object pointer
    #  @param idepic String depiction to be converted, str
    #  @param itype type of depiction provided as input, str
    #  @param otype types of depiction to be generated, {"", "", ..}
    #  @return odepic generated depictions, {"otype1": "odepic1", ..}
    def _convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise self.DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
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
        return odepic


    ## Function to create a dictionnary of old to new chemical id's
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkCIDdeprecated(self, cid):
        try:
            return self.deprecatedCID_cid[cid]
        except KeyError:
            return cid


    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkRIDdeprecated(self, rid):
        try:
            return self.deprecatedRID_rid[rid]
        except KeyError:
            return rid


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
    def _chebiXref(self):
        for cid in self.cid_xref:
            if 'chebi' in self.cid_xref[cid]:
                for c in self.cid_xref[cid]['chebi']:
                    self.chebi_cid[c] = cid


    def _inchikeyCID(self):
        for cid in self.cid_strc:
            if not self.cid_strc[cid]['inchikey'] in self.inchikey_cid:
                self.inchikey_cid[self.cid_strc[cid]['inchikey']] = []
            if not cid in self.inchikey_cid[self.cid_strc[cid]['inchikey']]:
                self.inchikey_cid[self.cid_strc[cid]['inchikey']].append(cid)


    #################################################################
    ################## Public functions #############################
    #################################################################


    ################## RetroRules parsers ####################################


    ## Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules. Uses the InchIkey as key to the dictionnary
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    #
    #  @param self Object pointer
    #  @param chem_prop_path Input file path
    def retroRulesStrc(self, rr_compounds_path):
        #TODO generate the structure when you get a chance
        for row in csv.DictReader(open(rr_compounds_path), delimiter='\t'):
            tmp = {'formula':  None,
                    'smiles': None,
                    'inchi': row['inchi'],
                    'inchikey': None}
            try:
                resConv = self._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles','inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except self.DepictionError as e:
                self.logger.warning('Could not convert some of the structures: '+str(tmp))
                self.logger.warning(e)
            if not self._checkCIDdeprecated(row['cid']) in self.cid_strc:
                self.cid_strc[self._checkCIDdeprecated(row['cid'])] = tmp


    ## Function to parse the rules_rall.tsv from RetroRules
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #
    #  @param self The object pointer.
    #  @param path The input file path.
    #  @return rule Dictionnary describing each reaction rule
    def retroReactions(self, rules_rall_path):
        try:
            #with open(rules_rall_path, 'r') as f:
            #    reader = csv.reader(f, delimiter = '\t')
            #    next(reader)
            #    rule = {}
            #    for row in reader:
            for row in csv.DictReader(open(rules_rall_path), delimiter='\t'):
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
                    self.rr_reactions[row['# Rule_ID']][row['Reaction_ID']] = {'rule_id': row['# Rule_ID'], 'rule_score': float(row['Score_normalized']), 'reac_id': self._checkRIDdeprecated(row['Reaction_ID']), 'subs_id': self._checkCIDdeprecated(row['Substrate_ID']), 'rel_direction': int(row['Rule_relative_direction']), 'left': {self._checkCIDdeprecated(row['Substrate_ID']): 1}, 'right': products}
                except ValueError:
                    self.logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    self.logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
        except FileNotFoundError as e:
                self.logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
                return {}


    ## Generate complete reactions from the rxn_recipes.tsv from RetroRules
    #
    #  These are the compplete reactions from which the reaction rules are generated from. This is used to
    # reconstruct the full reactions from monocomponent reactions
    #
    #  @param self The pointer object
    #  @param rxn_recipes_path Path to the recipes file
    #  @return Boolean that determines the success or failure of the function
    def retroRulesFullReac(self, rxn_recipes_path):
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


    ## Function to parse the comp_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of compartments id's (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    def mnxCompXref(self, comp_xref_path):
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
                        if not dbCompId in self.xref_comp:
                            self.xref_comp[dbCompId] = mnxc
        except FileNotFoundError:
            self.logger.error('comp_xref file not found')
            return {}



    #[TODO] merge the two functions
    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return Dictionnary of identifiers
    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
    def deprecatedMNXM(self, chem_xref_path):
        self.deprecatedCID_cid = {}
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


    ## Function to parse the reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param reac_xref_path Input file path
    #  @return Dictionnary of identifiers
    def deprecatedMNXR(self, reac_xref_path):
        self.deprecatedRID_rid = {}
        with open(reac_xref_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        self.deprecatedRID_rid[mnx[1]] = row[1]


    ## Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules
    # TODO: use the inchikey for all the cid
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    #
    #  @param self Object pointer
    #  @param chem_prop_path Input file path
    #  @return cid_strc Dictionnary of formula, smiles, inchi and inchikey
    def MNXstrc(self, chem_prop_path):
        with open(chem_prop_path) as f:
            c = csv.reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm = self._checkCIDdeprecated(row[0])
                    tmp = {'formula':  row[2],
                           'smiles': row[6],
                           'inchi': row[5],
                           'inchikey': row[8]}
                    for i in tmp:
                        if tmp[i]=='' or tmp[i]=='NA':
                            tmp[i] = None
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
                        try:
                            resConv = self._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                        except self.DepictionError as e:
                            self.logger.warning('Could not convert some of the structures: '+str(tmp))
                            self.logger.warning(e)
                        self.cid_strc[mnxm] = tmp


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
    def mnxXref(self, chem_xref_path):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser('Pass different Rules files')
    parser.add_argument('-rr_compounds_path', type=str, default=None)
    parser.add_argument('-rr_rules_path', type=str, default=None)
    parser.add_argument('-rr_rxn_recipes_path', type=str, default=None)
    params = parser.parse_args()
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
    rpcache = rpCache()
