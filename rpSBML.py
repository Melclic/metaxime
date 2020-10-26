"""rpSBML
.. moduleauthor:: Melchior du Lac
"""


import libsbml
from hashlib import md5
import os
import logging
import copy
import pandas as pd
import numpy as np


"""
TODO

1) fonction qui retourne l’identifiant du target “target_id"
2) fonction qui retourne directement la liste des identifiants des précurseurs “precursor_ids”.

"""

logging.basicConfig(
    level=logging.DEBUG,
    #level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)



## @package RetroPath SBML writer
# Documentation for SBML representation of the different model
#
# To exchange between the different workflow nodes, the SBML (XML) format is used. This
# implies using the libSBML library to create the standard definitions of species, reactions, etc...
# Here we also define our own annotations that are used internally in that we call BRSYNTH nodes.
# The object holds an SBML object and a series of methods to write and access BRSYNTH related annotations

##################################################################
############################### rpSBML ###########################
##################################################################


class rpSBML:
    """This class uses the libSBML object and handles it by adding BRSynth annotation
    """
    def __init__(self, model_name=None, document=None, path=None):
        """Constructor for the rpSBML class

        Note that the user can pass either a document libSBML object or a path to a SBML file. If a path is passed it overwrite the passed document object.

        :param model_name: The Name of the model
        :param document: The libSBML document class (Default: None)
        :param path: The path of a SBML file (Default: None)

        :type model_name: str
        :type path: str
        :type document: libsbml.SBMLDocument

        .. document private functions
        .. automethod:: _computeMeanRulesScore
        .. automethod:: _searchKey
        .. automethod:: _dictRPpathway
        .. automethod:: _findUniqueRowColumn
        .. automethod:: _checklibSBML
        .. automethod:: _nameToSbmlId
        .. automethod:: _genMetaID
        .. automethod:: _compareXref
        .. automethod:: _defaultBothAnnot
        .. automethod:: _defaultBRSynthAnnot
        .. automethod:: _defaultMIRIAMAnnot
        """
        self.logger = logging.getLogger(__name__)
        #WARNING: change this to reflect the different debugging levels
        self.logger.debug('Started instance of rpSBML')
        #self.logger.setLevel(logging.INFO)
        if model_name:
            self.model_name = model_name
        else:
            self.model_name = 'not_defined'
        self.mean_rules_score = 0.0
        if document==None and path==None:
            self.model = None
            self.document = None
        #document takes priority
        elif document:
            self.path = None
            self.model = self.document.getModel()
            #enabling the extra packages if they do not exists when reading a model
            if not self.model.isPackageEnabled('groups'):
                self._checklibSBML(self.model.enablePackage(
                    'http://www.sbml.org/sbml/level3/version1/groups/version1',
                    'groups',
                    True),
                        'Enabling the GROUPS package')
                self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
            if not self.model.isPackageEnabled('fbc'):
                self._checklibSBML(self.model.enablePackage(
                    'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                    'fbc',
                    True),
                        'Enabling the FBC package')
                self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')
            if self._isRPsbml():
                self.mean_rules_score = self._computeMeanRulesScore()
        elif path:
            self.path = path
            self.readSBML(path)
            if self._isRPsbml():
                self.mean_rules_score = self._computeMeanRulesScore()
        self.miriam_header = {'compartment': {'mnx': 'metanetx.compartment/', 'bigg': 'bigg.compartment/', 'seed': 'seed/', 'name': 'name/'}, 'reaction': {'mnx': 'metanetx.reaction/', 'rhea': 'rhea/', 'reactome': 'reactome/', 'bigg': 'bigg.reaction/', 'sabiork': 'sabiork.reaction/', 'ec': 'ec-code/', 'biocyc': 'biocyc/', 'lipidmaps': 'lipidmaps/', 'uniprot': 'uniprot/'}, 'species': {'inchikey': 'inchikey/', 'pubchem': 'pubchem.compound/','mnx': 'metanetx.chemical/', 'chebi': 'chebi/CHEBI:', 'bigg': 'bigg.metabolite/', 'hmdb': 'hmdb/', 'kegg_c': 'kegg.compound/', 'kegg_d': 'kegg.drug/', 'biocyc': 'biocyc/META:', 'seed': 'seed.compound/', 'metacyc': 'metacyc.compound/', 'sabiork': 'sabiork.compound/', 'reactome': 'reactome/R-ALL-'}}
        self.header_miriam = {'compartment': {'metanetx.compartment': 'mnx', 'bigg.compartment': 'bigg', 'seed': 'seed', 'name': 'name'}, 'reaction': {'metanetx.reaction': 'mnx', 'rhea': 'rhea', 'reactome': 'reactome', 'bigg.reaction': 'bigg', 'sabiork.reaction': 'sabiork', 'ec-code': 'ec', 'biocyc': 'biocyc', 'lipidmaps': 'lipidmaps', 'uniprot': 'uniprot'}, 'species': {'inchikey': 'inchikey', 'pubchem.compound': 'pubchem', 'metanetx.chemical': 'mnx', 'chebi': 'chebi', 'bigg.metabolite': 'bigg', 'hmdb': 'hmdb', 'kegg.compound': 'kegg_c', 'kegg.drug': 'kegg_d', 'biocyc': 'biocyc', 'seed.compound': 'seed', 'metacyc.compound': 'metacyc', 'sabiork.compound': 'sabiork', 'reactome': 'reactome'}}


    def __eq__(self, rpsbml):
        return \
            sorted(self.readRPpathwayIDs()) == sorted(rpsbml.readRPpathwayIDs()) \
        and self._dictRPpathway(self)==self._dictRPpathway(rpsbml)


    def __lt__(self, rpsbml):
        return self.mean_rules_score<rpsbml.mean_rules_score


    def __gt__(self, rpsbml):
        return self.mean_rules_score>rpsbml.mean_rules_score


    #######################################################################
    ########################## STATIC FUNCTIONS ###########################
    #######################################################################


    '''
    #TODO: remove this to another section of the project where more global are called
    @staticmethod
    def mergeSBMLFiles(path_source, path_target, path_merge):
        """Public function that merges two SBML files together

        :param path_source: Path of the source SBML file
        :param path_target: Path of the target SBML file
        :param path_merge: Path of the output SBML file

        :type path_source: str
        :type path_target: str
        :type path_merge: str

        :return: Success or failure of the function
        :rtype: bool
        """
        if not os.path.exists(path_target):
            logging.error('Target SBML file is invalid: '+str(path_target))
            return False
        if not os.path.exists(path_source):
            logging.error('Target SBML file is invalid: '+str(path_source))
            return False
        source_rpsbml = rpSBML('source', path=path_source)
        target_rpsbml = rpSBML('target', path=path_target)
        source_rpsbml.mergeModels(target_rpsbml.model)
        target_rpsbml.writeSBML(path_merge)
        return True
    '''

    #######################################################################
    ############################# PRIVATE FUNCTIONS ####################### 
    #######################################################################


    def _isRPsbml(self, pathway_id='rp_pathway'):
        """Check if the model is a rp_pathway
        """
        try:
            self.getGroupsMembers(pathway_id)
            return True
        except AttributeError:
            return False



    ############################# EXACT COMPARE ########################### 


    def _computeMeanRulesScore(self, pathway_id='rp_pathway'):
        """Compute the mean rules score

        :param pathway_id: The Groups id of the heterologous pathway

        :type pathway_id: str

        :rtype: float
        :return: The mean score of all the reaction rules
        """
        score = 0
        num_rules = 0
        for member in self.getGroupsMembers(pathway_id):
            reaction = self.model.getReaction(member)
            score += float(reaction.getAnnotation().getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild('rule_score').getAttrValue('value'))
            num_rules += 1
        if num_rules==0:
            return score
        else:
            return score/num_rules


    def _searchKey(self, keys, indict):
        for key in keys:
            if key in indict:
                return key


    def _dictRPpathway(self, rpsbml, pathway_id='rp_pathway'):
        """Put species in a dictionnary for further comparison

        :param rpsbml: rpSBML object
        :param pathway_id: The id of the Groups heterologous pathway

        :type rpsbml: rpSBML
        :type pathway_id: str

        :rtype: dict
        :return: Simplified dictionary version of rpSBML
        """
        #rpsbml = rpsbml.document.getModel()
        # Get Reactions
        reactions = {}
        for member in rpsbml.getGroupsMembers(pathway_id):
            reaction = rpsbml.model.getReaction(member)
            reactions[reaction.getId()] = self.readBRSYNTHAnnotation(reaction.getAnnotation())
        # Get Species
        species = {}
        for member in rpsbml.readUniqueRPspecies(pathway_id):
            spe = rpsbml.model.getSpecies(member)
            species[spe.getId()] = self.readBRSYNTHAnnotation(spe.getAnnotation())
        # Pathways dict
        d_reactions = {}
        keys = ['inchikey', 'inchi', 'smiles']
        # Select Reactions already loaded (w/o Sink one then)
        for reaction_id in reactions:
            # id = reactions[reaction]['smiles']
            #id = reaction_id
            d_reactions[reaction_id] = {}
            # Fill the reactants in a dedicated dict
            d_reactants = {}
            for reactant in rpsbml.model.getReaction(reaction_id).getListOfReactants():# inchikey / inchi sinon miriam sinon IDs
                # Il faut enregistrer toutes les infos (inchi, smiles, id)
                key = self._searchKey(keys, species[reactant.getSpecies()])
                if key:
                    key = species[reactant.getSpecies()][key]
                else:
                    key = reactant.getSpecies()
                d_reactants[key] = reactant.getStoichiometry()
            # Put all reactants dicts in reactions dict for which smiles notations are the keys
            d_reactions[reaction_id]['Reactants'] = d_reactants
            # Fill the products in a dedicated dict
            d_products = {}
            for product in rpsbml.model.getReaction(reaction_id).getListOfProducts():
                key = self._searchKey(keys, species[product.getSpecies()])
                if key:
                    key = species[product.getSpecies()][key]
                else:
                    key = product.getSpecies()
                d_products[key] = product.getStoichiometry()
            # Put all products dicts in reactions dict for which smiles notations are the keys
            d_reactions[reaction_id]['Products'] = d_products
        return d_reactions


    ############################ MERGE ####################################


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


    ############# OTHERS ##############


    def _checklibSBML(self, value, message):
        """Private function that checks the libSBML calls.

        Check that the libSBML python calls do not return error INT and if so, display the error. Taken from: http://sbml.org/Software/libSBML/docs/python-api/create_simple_model_8py-example.html

        :param value: The libSBML command returned int
        :param message: The string that describes the call

        :type value: int
        :type message: str

        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :return: None
        :rtype: None
        """
        if value is None:
            raise AttributeError('LibSBML returned a null value trying to ' + message + '.')
        elif type(value) is int:
            if value==libsbml.LIBSBML_OPERATION_SUCCESS:
                return
            else:
                err_msg = 'Error encountered trying to ' + message + '.' \
                        + 'LibSBML returned error code ' + str(value) + ': "' \
                        + libsbml.OperationReturnValue_toString(value).strip() + '"'
                raise AttributeError(err_msg)
        else:
            #self.logger.debug(message)
            return None


    def _nameToSbmlId(self, name):
        """String to SBML id's

        Convert any String to one that is compatible with the SBML meta_id formatting requirements

        :param name: The input string

        :type name: str

        :return: SBML valid string
        :rtype: str
        """
        IdStream = []
        count = 0
        end = len(name)
        if '0' <= name[count] and name[count] <= '9':
            IdStream.append('_')
        for count in range(0, end):
            if (('0' <= name[count] and name[count] <= '9') or
                    ('a' <= name[count] and name[count] <= 'z') or
                    ('A' <= name[count] and name[count] <= 'Z')):
                IdStream.append(name[count])
            else:
                IdStream.append('_')
        Id = ''.join(IdStream)
        if Id[len(Id) - 1] != '_':
            return Id
        return Id[:-1]


    def _genMetaID(self, name):
        """String to hashed id

        Hash an input string and then pass it to _nameToSbmlId()

        :param name: Input string

        :type name: str

        :return: Hashed string id
        :rtype: str
        """
        return self._nameToSbmlId(md5(str(name).encode('utf-8')).hexdigest())


    def _compareXref(self, current, toadd):
        """Compare two dictionaries of lists that describe the cross-reference and return the difference

        :param current: The source cross-reference dictionary
        :param toadd: The target cross-reference dictionary

        :type current: dict
        :type toadd: dict

        :return: Difference between the two cross-reference dictionaries
        :rtype: dict
        """
        toadd = copy.deepcopy(toadd)
        for database_id in current:
            try:
                list_diff = [i for i in toadd[database_id] if i not in current[database_id]]
                if not list_diff:
                    toadd.pop(database_id)
                else:
                    toadd[database_id] = list_diff
            except KeyError:
                pass
        if toadd==None:
            return []
        return toadd


    def _defaultBothAnnot(self, meta_id):
        """Returns a default annotation string that include MIRIAM and BRSynth annotation

        :param meta_id: The meta ID to be added to the default annotation

        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(meta_id or '')+'''">
      <bqbiol:is>
        <rdf:Bag>
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
    <rdf:BRSynth rdf:about="#'''+str(meta_id or '')+'''">
      <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
      </brsynth:brsynth>
    </rdf:BRSynth>
  </rdf:RDF>
</annotation>'''


    def _defaultBRSynthAnnot(self, meta_id):
        """Returns BRSynth default annotation string

        :param meta_id: The meta ID to be added to the annotation string

        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:BRSynth rdf:about="#'''+str(meta_id or '')+'''">
      <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
      </brsynth:brsynth>
    </rdf:BRSynth>
  </rdf:RDF>
</annotation>'''


    def _defaultMIRIAMAnnot(self, meta_id):
        """Returns MIRIAM default annotation string

        :param meta_id: The meta ID to be added to the annotation string

        :type meta_id: str

        :return: The default annotation string
        :rtype: str
        """
        return '''<annotation>
  <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/">
    <rdf:Description rdf:about="#'''+str(meta_id or '')+'''">
      <bqbiol:is>
        <rdf:Bag>
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
  </rdf:RDF>
</annotation>'''


    def updateBRSynthPathway(self, rpsbml_dict, pathway_id='rp_pathway'):
        """Given a dict of the same as rpsbml_dict (output of rpSBML.asdict), update the differences

        :param rpsbml_dict: Heterologous dict
        :param pathway_id: The Groups id of the heterologous pathway

        :type rpsbml_dict: dict
        :type pathway_id: str

        :rtype: None
        :return: None
        """
        self.logger.debug('rpsbml_dict: '+str(rpsbml_dict))
        groups = self.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        for bd_id in rpsbml_dict['pathway']['brsynth']:
            if bd_id[:5]=='norm_' or bd_id=='global_score':
                try:
                    value = rpsbml_dict['pathway']['brsynth'][bd_id]['value']
                except KeyError:
                    self.logger.warning('The entry '+str(db_id)+' doesnt contain value')
                    self.logger.warning('No" value", using the root')
                try:
                    units = rpsbml_dict['pathway']['brsynth'][bd_id]['units']
                except KeyError:
                    units = None
                self.addUpdateBRSynth(rp_pathway, bd_id, value, units, False)
        for reac_id in rpsbml_dict['reactions']:
            reaction = self.model.getReaction(reac_id)
            if reaction==None:
                self.logger.warning('Skipping updating '+str(reac_id)+', cannot retreive it')
                continue
            for bd_id in rpsbml_dict['reactions'][reac_id]['brsynth']:
                if bd_id[:5]=='norm_':
                    try:
                        value = rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['value']
                    except KeyError:
                        value = rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]
                        self.logger.warning('No" value", using the root')
                    try:
                        units = rpsbml_dict['reactions'][reac_id]['brsynth'][bd_id]['units']
                    except KeyError:
                        units = None
                    self.addUpdateBRSynth(reaction, bd_id, value, units, False)


    ######################################################################
    ####################### PUBLIC FUNCTIONS #############################
    ######################################################################


    def addMIRIAMinchiKey(self):
        """Check the MIRIAM annotation for MetaNetX or CHEBI id's and try to recover the inchikey from cache and add it to MIRIAM

        :rtype: bool

        :return: Success or failure of the function
        """
        for spe in self.model.getListOfSpecies():
            inchikey = None
            miriam_dict = self.readMIRIAMAnnotation(spe.getAnnotation())
            if 'inchikey' in miriam_dict:
                self.logger.debug('The species '+str(spe.id)+' already has an inchikey... skipping')
                continue
            try:
                for mnx in miriam_dict['metanetx']:
                    inchikey = self.cid_strc[self._checkCIDdeprecated(mnx)]['inchikey']
                    if inchikey:
                        self.addUpdateMIRIAM(spe, 'species', {'inchikey': [inchikey]})
                    else:
                        self.logger.warning('The inchikey is empty for: '+str(spe.id))
                    continue
            except KeyError:
                try:
                    for chebi in miriam_dict['chebi']:
                        inchikey = self.cid_strc[self._checkCIDdeprecated(self.chebi_cid[chebi])]['inchikey']
                        if inchikey:
                            self.addUpdateMIRIAM(spe, 'species', {'inchikey': [inchikey]})
                        else:
                            self.logger.warning('The inchikey is empty for: '+str(spe.id))
                        continue
                except KeyError:
                    self.logger.warning('Cannot find the inchikey for: '+str(spe.id))
        return True


    ############################# MERGE FUNCTIONS #########################

    #REACTION

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

    #SPECIES

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

    #EC NUMBER

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

    #MERGE


    #TODO: add a confidence in the merge using the score in 
    #TODO: seperate the different parts so that others may use it
    def mergeModels(self, input_model, output_path=None):
        """Merge two models species and reactions using the annotations to recognise the same species and reactions

        The source model has to have both the GROUPS and FBC packages enabled in its SBML. The course must have a groups
        called rp_pathway. If not use the readSBML() function to create a model
        We add the reactions and species from the rpsbml to the target_model. When the output_path is specified, the function
        will output the model to that path and if left to the default None, the self model will be overwritten.

        :param target_model: The target model object
        :param output_path: The path to the output file. Note that if its set to None, then self model is overwritten

        :type target_model: libsbml.Model
        :type output_path: None

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
        targetParametersID = [i.getId() for i in target_model.getListOfParameters()]
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
        '''
        if fillOrphanSpecies==True:
            self.fillOrphan(target_rpsbml, self.pathway_id, compartment_id)
        '''
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

    
    #TODO: need to pass the species and reaction conversions that the above function returns
    def overwriteRPannot(self, source_rpsbml, species_source_target, reactions_source_target, source_pathway_id='rp_pathway', pathway_id='rp_pathway'):
        """Given a rpsbml of entry, overwrite the current annotations of a pathway. This is useful when performing operations in merged model, i.e. rpFBA and then
        adding this info to the new model

        :param source_rpsbml: The source rpSBML from which to read the annotations of the pathway
        :param source_pathway_id: The Groups id of the source heterologous pathway
        :param pathway_id: The Groups id of the heterologous pathway

        :type source_rpsbml: rpSBML
        :param source_pathway_id: str 
        :param pathway_id: str 

        :rtype: None
        :return: None
        """
        self.logger.debug('Returning model with heterologous pathway only')
        rev_reactions_convert = {v: k for k, v in reactions_convert.items()}
        #### species #####
        #TODO: overwrite all the species:
        groups = source_rpsbml.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        self.logger.debug('---- Reactions ----')
        for member in rp_pathway.getListOfMembers():
            #### reaction annotation
            self.logger.debug(member.getIdRef())
            reac_fba = source_rpsbml.model.getReaction(member.getIdRef())
            self.logger.debug(reac_fba)
            try:
                #reac_in = self.model.getReaction(reactions_convert[member.getIdRef()])
                reac_in = self.model.getReaction(rev_reactions_convert[member.getIdRef()])
            except KeyError:
                reac_in = self.model.getReaction(member.getIdRef())
            self.logger.debug(reac_in)
            self.logger.debug(reac_fba.getAnnotation())
            reac_in.setAnnotation(reac_fba.getAnnotation())
            #### species TODO: only for shadow price
        #### add groups ####
        source_groups = source_rpsbml.model.getPlugin('groups')
        target_groups = self.model.getPlugin('groups')
        target_groups_id = [i.getId() for i in target_groups.getListOfGroups()]
        for source_group in source_groups.getListOfGroups():
            #self.logger.debug('Replacing group id: '+str(source_group.getId()))
            if source_group.getId()==species_group_id:
                target_group = target_groups.getGroup(source_group.getId())
                #TODO: #### replace the new potentially incorect central species with the normal ones #####
                #delete all the previous members
                self.logger.debug('Removing central_species')
                for i in range(target_group.getNumMembers()):
                    self.logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                #add the new ones
                for cs in cent_spe:
                    self.logger.debug('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)
            elif source_group.getId()==sink_species_group_id:
                target_group = target_groups.getGroup(source_group.getId())
                self.logger.debug('Removing sink species')
                for i in range(target_group.getNumMembers()):
                    self.logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                #add the new ones
                for cs in sink_spe:
                    self.logger.debug('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)
            elif source_group.getId() in target_groups_id:
                target_group = target_groups.getGroup(source_group.getId())
                target_group.setAnnotation(source_group.getAnnotation())
            '''
            elif source_group.getId()==pathway_id:
                target_group = target_groups.getGroup(source_group.getId())
                self.logger.debug('Removing rp ractions')
                for i in range(target_group.getNumMembers()):
                    self.logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                #add the new ones
                for cs in rp_reac:
                    self.logger.debug('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)  
            '''
        #### add objectives ####
        source_fbc = source_rpsbml.model.getPlugin('fbc')
        target_fbc = self.model.getPlugin('fbc')
        target_objID = [i.getId() for i in target_fbc.getListOfObjectives()]
        for source_obj in source_fbc.getListOfObjectives():
            source_obj_id = source_obj.getId()
            if source_obj.getId() in target_objID:
                target_obj = target_fbc.getObjective(source_obj.getId())
                target_obj.setAnnotation(source_obj.getAnnotation())
                for target_fluxObj in target_obj.getListOfFluxObjectives():
                    for source_fluxObj in source_obj.getListOfFluxObjectives():
                        if target_fluxObj.getReaction()==source_fluxObj.getReaction():
                            target_fluxObj.setAnnotation(source_fluxObj.getAnnotation())
            else:
                target_fbc.addObjective(source_obj)
        #rpsbml.createMultiFluxObj('obj_RP1_sink', ['RP1_sink'], [1])
        target_fbc.setActiveObjectiveId(source_obj_id) #tmp random assigenement of objective
 


    ####################### ADD ELEMENTS ############


    def addUpdateBRSynth(self, sbase_obj, annot_header, value, units=None, isAlone=False, isList=False, isSort=True, meta_id=None):
        """Append or update an entry to the BRSynth annotation of the passed libsbml.SBase object.

        If the annot_header isn't contained in the annotation it is created. If it already exists it overwrites it

        :param sbase_obj: The libSBML object to add the different
        :param annot_header: The annotation header that defines the type of entry
        :param value: The value(s) to add
        :param units: Add a values unit to the entry
        :param isAlone: Add the entry without any units or defined within a value child (Setting this to True will ignore any units)
        :param isList: Define if the value entry is a list or not
        :param isSort: Sort the list that is passed (Only if the isList is True)
        :param meta_id: The meta ID to be added to the annotation string

        :type sbase_obj: libsbml.SBase
        :type annot_header: str
        :type value: Union[str, int, float, list]
        :type units: str
        :type isAlone: bool
        :type isList: bool
        :type isSort: bool
        :type meta_id: str

        :rtype: bool
        :return: Sucess or failure of the function
        """
        self.logger.debug('############### '+str(annot_header)+' ################')
        if isList:
            annotation = '''<annotation>
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
        <rdf:BRSynth rdf:about="#adding">
          <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">
            <brsynth:'''+str(annot_header)+'''>'''
            if isSort:
                for name in sorted(value, key=value.get, reverse=True):
                    if isAlone:
                        annotation += '<brsynth:'+str(name)+'>'+str(value[name])+'</brsynth:'+str(name)+'>'
                    else:
                        if units:
                            annotation += '<brsynth:'+str(name)+' units="'+str(units)+'" value="'+str(value[name])+'" />'
                        else:
                            annotation += '<brsynth:'+str(name)+' value="'+str(value[name])+'" />'
            else:
                for name in value:
                    if isAlone:
                        annotation += '<brsynth:'+str(name)+'>'+str(value[name])+'</brsynth:'+str(name)+'>'
                    else:
                        if units:
                            annotation += '<brsynth:'+str(name)+' units="'+str(units)+'" value="'+str(value[name])+'" />'
                        else:
                            annotation += '<brsynth:'+str(name)+' value="'+str(value[name])+'" />'
            annotation += '''
            </brsynth:'''+str(annot_header)+'''>
          </brsynth:brsynth>
        </rdf:BRSynth>
      </rdf:RDF>
    </annotation>'''
        else:
            #### create the string
            annotation = '''<annotation>
      <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
        <rdf:BRSynth rdf:about="#adding">
          <brsynth:brsynth xmlns:brsynth="http://brsynth.eu">'''
            if isAlone:
                annotation += '<brsynth:'+str(annot_header)+'>'+str(value)+'</brsynth:'+str(annot_header)+'>'
            else:
                if units:
                    annotation += '<brsynth:'+str(annot_header)+' units="'+str(units)+'" value="'+str(value)+'" />'
                else:
                    annotation += '<brsynth:'+str(annot_header)+' value="'+str(value)+'" />'
            annotation += '''
          </brsynth:brsynth>
        </rdf:BRSynth>
      </rdf:RDF>
    </annotation>'''
        annot_obj = libsbml.XMLNode.convertStringToXMLNode(annotation)
        if annot_obj==None:
            self.logger.error('Cannot conver this string to annotation object: '+str(annotation))
            return False
        #### retreive the annotation object
        brsynth_annot = None
        obj_annot = sbase_obj.getAnnotation()
        if not obj_annot:
            sbase_obj.setAnnotation(libsbml.XMLNode.convertStringToXMLNode(self._defaultBRSynthAnnot(meta_id)))
            obj_annot = sbase_obj.getAnnotation()
            if not obj_annot:
                self.logger.error('Cannot update BRSynth annotation')
                return False
        brsynth_annot = obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        if not brsynth_annot:
             self.logger.error('Cannot find the BRSynth annotation')
             return False
        #add the annotation and replace if it exists
        isfound_target = False
        #self.logger.debug(brsynth_annot.toXMLString())
        for i in range(brsynth_annot.getNumChildren()):
            self.logger.debug(annot_header+' -- '+str(brsynth_annot.getChild(i).getName()))
            if annot_header==brsynth_annot.getChild(i).getName():
                isfound_target = True
                '''
                self._checklibSBML(brsynth_annot.removeChild(brsynth_annot.getIndex(i)),
                    'Removing annotation '+str(annot_header))
                '''
                self._checklibSBML(brsynth_annot.removeChild(i), 'Removing annotation '+str(annot_header))
                isfound_source = False
                source_brsynth_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth')
                for y in range(source_brsynth_annot.getNumChildren()):
                    self.logger.debug('\t'+annot_header+' -- '+str(source_brsynth_annot.getChild(y).getName()))
                    if str(annot_header)==str(source_brsynth_annot.getChild(y).getName()):
                        isfound_source = True
                        self.logger.debug('Adding annotation to the brsynth annotation: '+str(source_brsynth_annot.getChild(y).toXMLString()))
                        towrite_annot = source_brsynth_annot.getChild(y)
                        self.logger.debug(brsynth_annot.toXMLString())
                        self._checklibSBML(brsynth_annot.addChild(towrite_annot), ' 1 - Adding annotation to the brsynth annotation')
                        self.logger.debug(brsynth_annot.toXMLString())
                        break
                if not isfound_source:
                    self.logger.error('Cannot find '+str(annot_header)+' in source annotation')
        if not isfound_target:
            self.logger.debug('Cannot find '+str(annot_header)+' in target annotation')
            isfound_source = False
            source_brsynth_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth')
            for y in range(source_brsynth_annot.getNumChildren()):
                self.logger.debug('\t'+annot_header+' -- '+str(source_brsynth_annot.getChild(y).getName()))
                if str(annot_header)==str(source_brsynth_annot.getChild(y).getName()):
                    isfound_source = True
                    self.logger.debug('Adding annotation to the brsynth annotation: '+str(source_brsynth_annot.getChild(y).toXMLString()))
                    towrite_annot = source_brsynth_annot.getChild(y)
                    self.logger.debug(brsynth_annot.toXMLString())
                    self._checklibSBML(brsynth_annot.addChild(towrite_annot), '2 - Adding annotation to the brsynth annotation')
                    self.logger.debug(brsynth_annot.toXMLString())
                    break
            if not isfound_source:
                self.logger.error('Cannot find '+str(annot_header)+' in source annotation')
            #toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
            #self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
                return False
        '''
        if brsynth_annot.getChild(annot_header).toXMLString()=='':
            toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
            self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
        else:
            #try:
            self.logger.debug('==============================')
            found_child = False
            for i in range(brsynth_annot.getNumChildren()):
                if annot_header==brsynth_annot.getChild(i).getName():
                    self.logger.debug('Found the same name to remove: '+str(annot_header))
                    self._checklibSBML(brsynth_annot.removeChild(brsynth_annot.getIndex(i)),
                        'Removing annotation '+str(annot_header))
                    toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
                    self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
                    found_child = True
                    break
            #cause by a bbug with string lookup
            if not found_child:
                self.logger.warning('Bug with lookup adding it now: '+str(annot_header))
                toWrite_annot = annot_obj.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(annot_header)
                self._checklibSBML(brsynth_annot.addChild(toWrite_annot), 'Adding annotation to the brsynth annotation')
            #except OverflowError:
            #    self.logger.warning('TODO: Overflow error that must be dealt with')
            #    self.logger.warning(brsynth_annot.getChild(annot_header).toXMLString())
            #    return False
        '''
        return True


    def addUpdateMIRIAM(self, sbase_obj, type_param, xref, meta_id=None):
        """Append or update an entry to the MIRIAM annotation of the passed libsbml.SBase object.

        If the annot_header isn't contained in the annotation it is created. If it already exists it overwrites it

        :param sbase_obj: The libSBML object to add the different 
        :param type_param: The type of parameter entered. Valid include ['compartment', 'reaction', 'species']
        :param xref: Dictionnary of the cross reference
        :param meta_id: The meta ID to be added to the annotation string

        :type sbase_obj: libsbml.SBase
        :type type_param: str
        :type xref: dict
        :type meta_id: str

        :rtype: bool
        :return: Sucess or failure of the function
        """
        if not type_param in ['compartment', 'reaction', 'species']:
            self.logger.error('type_param must be '+str(['compartment', 'reaction', 'species'])+' not '+str(type_param))
            return False
        miriam_annot = None
        isReplace = False
        try:
            miriam_annot = sbase_obj.getAnnotation().getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            miriam_elements = self.readMIRIAMAnnotation(sbase_obj.getAnnotation())
            if not miriam_elements:
                isReplace = True
                if not meta_id:
                    meta_id = self._genMetaID('tmp_addUpdateMIRIAM')
                miriam_annot_1 = libsbml.XMLNode.convertStringToXMLNode(self._defaultBothAnnot(meta_id))
                miriam_annot = miriam_annot_1.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            else:
                miriam_elements = None
        except AttributeError:
            try:
                #Cannot find MIRIAM annotation, create it
                isReplace = True
                if not meta_id:
                    meta_id = self._genMetaID('tmp_addUpdateMIRIAM')
                miriam_annot = libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id))
                miriam_annot = miriam_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            except AttributeError:
                self.logger.error('Fatal error fetching the annotation')
                return False
        #compile the list of current species
        inside = {}
        for i in range(miriam_annot.getNumChildren()):
            single_miriam = miriam_annot.getChild(i)
            if single_miriam.getAttributes().getLength()>1:
                self.logger.error('MIRIAM annotations should never have more than 1: '+str(single_miriam.toXMLString()))
                continue
            single_miriam_attr = single_miriam.getAttributes()
            if not single_miriam_attr.isEmpty():
                try:
                    db = single_miriam_attr.getValue(0).split('/')[-2]
                    v = single_miriam_attr.getValue(0).split('/')[-1]
                    inside[self.header_miriam[type_param][db]].append(v)
                except KeyError:
                    try:
                        db = single_miriam_attr.getValue(0).split('/')[-2]
                        v = single_miriam_attr.getValue(0).split('/')[-1]
                        inside[self.header_miriam[type_param][db]] = [v]
                    except KeyError:
                        self.logger.warning('Cannot find the self.header_miriram entry '+str(db))
                        continue
            else:
                self.logger.warning('Cannot return MIRIAM attribute')
                pass
        #add or ignore
        self.logger.debug('inside: '+str(inside))
        self.logger.debug('xref: '+str(xref))
        toadd = self._compareXref(inside, xref)
        self.logger.debug('toadd: '+str(toadd))
        for database_id in toadd:
            for species_id in toadd[database_id]:
                #not sure how to avoid having it that way
                if database_id in self.miriam_header[type_param]:
                    try:
                        #determine if the dictionnaries
                        annotation = '''<annotation>
    <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
    <rdf:Description rdf:about="#tmp">
      <bqbiol:is>
        <rdf:Bag>'''
                        if type_param=='species':
                            if database_id=='kegg' and species_id[0]=='C':
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param]['kegg_c']+str(species_id)+'''"/>'''
                            elif database_id=='kegg' and species_id[0]=='D':
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param]['kegg_d']+str(species_id)+'''"/>'''
                            else:
                                annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param][database_id]+str(species_id)+'''"/>'''
                        else:
                            annotation += '''
              <rdf:li rdf:resource="http://identifiers.org/'''+self.miriam_header[type_param][database_id]+str(species_id)+'''"/>'''
                        annotation += '''
        </rdf:Bag>
      </bqbiol:is>
    </rdf:Description>
    </rdf:RDF>
    </annotation>'''
                        toPass_annot = libsbml.XMLNode.convertStringToXMLNode(annotation)
                        toWrite_annot = toPass_annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag').getChild(0)
                        miriam_annot.insertChild(0, toWrite_annot)
                    except KeyError:
                        #WARNING need to check this
                        self.logger.warning('Cannot find '+str(database_id)+' in self.miriam_header for '+str(type_param))
                        continue
        if isReplace:
            ori_miriam_annot = sbase_obj.getAnnotation()
            if ori_miriam_annot==None:
                sbase_obj.unsetAnnotation()
                sbase_obj.setAnnotation(miriam_annot)
            else:
                self._checklibSBML(ori_miriam_annot.getChild('RDF').getChild('Description').getChild('is').removeChild(0), 'Removing annotation "is"')
                self._checklibSBML(ori_miriam_annot.getChild('RDF').getChild('Description').getChild('is').addChild(miriam_annot), 'Adding annotation to the brsynth annotation')
        return True


    def asdict(self, pathway_id='rp_pathway'):
        """Generate the dictionnary of all the annotations of a pathway species, reaction and pathway annotations

        :param pathway_id: The pathway ID (Default: rp_pathway)

        :type pathway_id: str

        :rtype: dict
        :return: Dictionnary of the pathway annotation
        """
        groups = self.model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        reactions = rp_pathway.getListOfMembers()
        #pathway
        rpsbml_json = {}
        rpsbml_json['pathway'] = {}
        rpsbml_json['pathway']['brsynth'] = self.readBRSYNTHAnnotation(rp_pathway.getAnnotation())
        #reactions
        rpsbml_json['reactions'] = {}
        for member in reactions:
            reaction = self.model.getReaction(member.getIdRef())
            annot = reaction.getAnnotation()
            rpsbml_json['reactions'][member.getIdRef()] = {}
            rpsbml_json['reactions'][member.getIdRef()]['brsynth'] = self.readBRSYNTHAnnotation(annot)
            rpsbml_json['reactions'][member.getIdRef()]['miriam'] = self.readMIRIAMAnnotation(annot)
        #loop though all the species
        rpsbml_json['species'] = {}
        for spe_id in self.readUniqueRPspecies(pathway_id):
            species = self.model.getSpecies(spe_id)
            annot = species.getAnnotation()
            rpsbml_json['species'][spe_id] = {}
            rpsbml_json['species'][spe_id]['brsynth'] = self.readBRSYNTHAnnotation(annot)
            rpsbml_json['species'][spe_id]['miriam'] = self.readMIRIAMAnnotation(annot)
        return rpsbml_json


    ########################## INPUT/OUTPUT #############################


    def readSBML(self, path, model_name=None):
        """Open an SBML file to the object

        :param path: Path to the input SBML file

        :type path: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: None
        :return: Dictionnary of the pathway annotation
        """
        if not os.path.isfile(path):
            self.logger.error('Invalid input file')
            raise FileNotFoundError
        document = libsbml.readSBMLFromFile(path)
        self._checklibSBML(document, 'reading input file')
        errors = document.getNumErrors()
        #display the errors in the log accordning to the severity
        for err in [document.getError(i) for i in range(document.getNumErrors())]:
            #TODO if the error is related to packages not enabled (like groups or fbc) activate them
            if err.isFatal:
                self.logger.error('libSBML reading error: '+str(err.getShortMessage()))
                raise FileNotFoundError
            else:
                self.logger.warning('libSBML reading warning: '+str(err.getShortMessage()))
        model = document.getModel()
        if not model:
            self.logger.error('Either the file was not read correctly or the SBML is empty')
            raise FileNotFoundError
        self.document = document
        self.model = model
        if self._isRPsbml():
            self.mean_rules_score = self._computeMeanRulesScore()
        if model_name:
            self.model_name = model_name
        elif not self.model_name:
            self.model_name = model.getId().replace('model_', '').lower()
        #enabling the extra packages if they do not exists when reading a model
        if not self.model.isPackageEnabled('groups'):
            self._checklibSBML(self.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/groups/version1',
                'groups',
                True),
                    'Enabling the GROUPS package')
            self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        if not self.model.isPackageEnabled('fbc'):
            self._checklibSBML(self.model.enablePackage(
                'http://www.sbml.org/sbml/level3/version1/fbc/version2',
                'fbc',
                True),
                    'Enabling the FBC package')
            self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')


    def writeSBML(self, path=None):
        """Export the metabolic network to a SBML file

        :param path: Path to the output SBML file

        :type path: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: bool
        :return: Success or failure of the command
        """
        ####### check the path #########
        #need to determine where are the path id's coming from
        try:
            p = None
            if path:
                if os.path.isdir(path):
                    p = os.path.join(path, str(self.model_name)+'_rpsbml.xml')
                else:
                    if os.path.exists(path):
                        self.logger.warning('overwriting the following sbml file: '+str(path))
                    p = path
            else:
                if self.path:
                    self.logger.warning('overwriting the following sbml file: '+str(self.path))
                    p = self.path
                else:
                    self.logger.error('Need to define a file or folder output')
                    return False
            ########## check and create folder #####
        except FileNotFoundError:
            self.logger.error('Cannot find the following path: '+str(path))
            return False
        libsbml.writeSBMLToFile(self.document, p)
        return True


    ########################## FindCreate ###############################


    #TODO: seperate the find and the create
    def findCreateObjective(self, reactions, coefficients, is_max=True, objective_id=None):
        """Find the objective (with only one reaction associated) based on the reaction ID and if not found create it

        :param reactions: List of the reactions id's to set as objectives
        :param coefficients: List of the coefficients about the objectives
        :param is_max: Maximise or minimise the objective
        :param objective_id: overwite the default id if created (from obj_[reactions])

        :type reactions: list
        :type coefficients: list
        :type is_max: bool
        :type objective_id: str

        :raises FileNotFoundError: If the file cannot be found
        :raises AttributeError: If the libSBML command encounters an error or the input value is None

        :rtype: str
        :return: Objective ID
        """
        fbc_plugin = self.model.getPlugin('fbc')
        self._checklibSBML(fbc_plugin, 'Getting FBC package')
        if not objective_id:
            objective_id = 'obj_'+'_'.join(reactions)
            self.logger.debug('Setting objective as '+str(objective_id))
        for objective in fbc_plugin.getListOfObjectives():
            if objective.getId()==objective_id:
                self.logger.warning('The specified objective id ('+str(objective_id)+') already exists')
                return objective_id
            if not set([i.getReaction() for i in objective.getListOfFluxObjectives()])-set(reactions):
                #TODO: consider setting changing the name of the objective
                self.logger.warning('The specified objective id ('+str(objective_id)+') has another objective with the same reactions: '+str(objective.getId()))
                return objective.getId()
        #If cannot find a valid objective create it
        self.createMultiFluxObj(objective_id,
                                reactions,
                                coefficients,
                                is_max)
        return objective_id


    ########################## READ #####################################


    def readRPpathwayIDs(self, pathway_id='rp_pathway'):
        """Return the members of rp_pathways

        NOTE: Try not to use this function and instead use getGroupsMembers()

        :rtype: list
        :return: List of member id's of a particular group
        """
        groups = self.model.getPlugin('groups')
        groups_obj = groups.getGroup(pathway_id)
        self._checklibSBML(groups_obj, 'retreiving groups groups_obj (readRPpathwayIDs): '+str(group_id))
        toRet = []
        for member in groups_obj.getListOfMembers():
            toRet.append(member.getIdRef())
        return toRet


    #TODO: rename this function to readGroupsMembers
    #TODO: add error handling if the groups does not exist
    #TODO: change the pathway_id to groups_id
    def getGroupsMembers(self, group_id):
        """Return the members of a groups entry

        :param group_id: The pathway ID

        :type group_id: str

        :rtype: list
        :return: List of member id's of a particular group
        """
        groups = self.model.getPlugin('groups')
        groups_obj = groups.getGroup(group_id)
        self._checklibSBML(groups_obj, 'retreiving groups groups_obj (getGroupsMembers): '+str(group_id))
        toRet = []
        for member in groups_obj.getListOfMembers():
            toRet.append(member.getIdRef())
        return toRet


    def readRPrules(self, pathway_id='rp_pathway'):
        """Return the list of reaction rules contained within a pathway

        :param pathway_id: The pathway ID (Default: rp_pathway)

        :type pathway_id: str

        :rtype: dict
        :return: Dictionnary of reaction rules (rule_id as key)
        """
        toRet = {}
        for reacId in self.getGroupsMembers(pathway_id):
            reac = self.model.getReaction(reacId)
            brsynth_annot = self.readBRSYNTHAnnotation(reac.getAnnotation())
            if not brsynth_annot['rule_id']=='' and not brsynth_annot['smiles']=='':
                toRet[brsynth_annot['rule_id']] = brsynth_annot['smiles'].replace('&gt;', '>')
        return toRet


    #TODO: merge with unique species
    #TODO: change the name of the function to read
    def readRPspecies(self, pathway_id='rp_pathway'):
        """Return the species stoichiometry of a pathway

        :param pathway_id: The pathway ID (Default: rp_pathway)

        :type pathway_id: str

        :rtype: dict
        :return: Dictionary of the pathway species and reactions
        """
        reacMembers = {}
        for reacId in self.getGroupsMembers(pathway_id):
            reacMembers[reacId] = {}
            reacMembers[reacId]['products'] = {}
            reacMembers[reacId]['reactants'] = {}
            reac = self.model.getReaction(reacId)
            for pro in reac.getListOfProducts():
                reacMembers[reacId]['products'][pro.getSpecies()] = pro.getStoichiometry()
            for rea in reac.getListOfReactants():
                reacMembers[reacId]['reactants'][rea.getSpecies()] = rea.getStoichiometry()
        return reacMembers


    def readUniqueRPspecies(self, pathway_id='rp_pathway'):
        """Return the unique species of a pathway

        :param pathway_id: The pathway ID (Default: rp_pathway)

        :type pathway_id: str

        :rtype: list
        :return: List of unique species
        """
        rpSpecies = self.readRPspecies()
        toRet = []
        for i in rpSpecies:
            for y in rpSpecies[i]:
                for z in rpSpecies[i][y]:
                    if not z in toRet:
                        toRet.append(z)
        return toRet
        #reacMembers = self.readRPspecies(pathway_id)
        #return set(set(ori_rp_path['products'].keys())|set(ori_rp_path['reactants'].keys()))


    def readTaxonAnnotation(self, annot):
        """Return he taxonomy ID from an annotation

        :param annot: The annotation object of libSBML

        :type annot: libsbml.XMLNode

        :rtype: dict
        :return: Dictionary of all taxonomy id's
        """
        try:
            toRet = {}
            bag = annot.getChild('RDF').getChild('Description').getChild('hasTaxon').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                if str_annot=='':
                    self.logger.warning('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                    continue
                dbid = str_annot.split('/')[-2].split('.')[0]
                if len(str_annot.split('/')[-1].split(':'))==2:
                    cid = str_annot.split('/')[-1].split(':')[1]
                else:
                    cid = str_annot.split('/')[-1]
                if not dbid in toRet:
                    toRet[dbid] = []
                toRet[dbid].append(cid)
            return toRet
        except AttributeError:
            return {}


    def readMIRIAMAnnotation(self, annot):
        """Return the MIRIAM annotations of species

        :param annot: The annotation object of libSBML

        :type annot: libsbml.XMLNode

        :rtype: dict
        :return: Dictionary of all the annotation of species
        """
        try:
            toRet = {}
            bag = annot.getChild('RDF').getChild('Description').getChild('is').getChild('Bag')
            for i in range(bag.getNumChildren()):
                str_annot = bag.getChild(i).getAttrValue(0)
                if str_annot=='':
                    self.logger.warning('This contains no attributes: '+str(bag.getChild(i).toXMLString()))
                    continue
                dbid = str_annot.split('/')[-2].split('.')[0]
                if len(str_annot.split('/')[-1].split(':'))==2:
                    cid = str_annot.split('/')[-1].split(':')[1]
                else:
                    cid = str_annot.split('/')[-1]
                if not dbid in toRet:
                    toRet[dbid] = []
                toRet[dbid].append(cid)
            return toRet
        except AttributeError:
            return {}


    def readBRSYNTHAnnotation(self, annot):
        """Return a dictionnary of all the information in a BRSynth annotations

        :param annot: The annotation object of libSBML

        :type annot: libsbml.XMLNode

        :rtype: dict
        :return: Dictionary of all the BRSynth annotations
        """
        toRet = {'dfG_prime_m': None,
                 'dfG_uncert': None,
                 'dfG_prime_o': None,
                 'path_id': None,
                 'step_id': None,
                 'sub_step_id': None,
                 'rule_score': None,
                 'smiles': None,
                 'inchi': None,
                 'inchikey': None,
                 'selenzyme': None,
                 'rule_id': None,
                 'rule_ori_reac': None,
                 'rule_score': None,
                 'global_score': None}
        if annot==None:
            self.logger.warning('The passed annotation is None')
            return {}
        bag = annot.getChild('RDF').getChild('BRSynth').getChild('brsynth')
        for i in range(bag.getNumChildren()):
            ann = bag.getChild(i)
            if ann=='':
                self.logger.warning('This contains no attributes: '+str(ann.toXMLString()))
                continue
            if ann.getName()=='dfG_prime_m' or ann.getName()=='dfG_uncert' or ann.getName()=='dfG_prime_o' or ann.getName()[0:4]=='fba_' or ann.getName()=='flux_value':
                try:
                    toRet[ann.getName()] = {
                            'units': ann.getAttrValue('units'),
                            'value': float(ann.getAttrValue('value'))}
                except ValueError:
                    self.logger.warning('Cannot interpret '+str(ann.getName())+': '+str(ann.getAttrValue('value')+' - '+str(ann.getAttrValue('units'))))
                    toRet[ann.getName()] = {
                            'units': None,
                            'value': None}
            elif ann.getName()=='path_id' or ann.getName()=='step_id' or ann.getName()=='sub_step_id':
                try:
                    #toRet[ann.getName()] = int(ann.getAttrValue('value'))
                    toRet[ann.getName()] = {'value': int(ann.getAttrValue('value'))}
                except ValueError:
                    toRet[ann.getName()] = None
            elif ann.getName()=='rule_score' or ann.getName()=='global_score' or ann.getName()[:5]=='norm_':
                try:
                    #toRet[ann.getName()] = float(ann.getAttrValue('value'))
                    toRet[ann.getName()] = {'value': float(ann.getAttrValue('value'))}
                except ValueError:
                    toRet[ann.getName()] = None
            elif ann.getName()=='smiles':
                toRet[ann.getName()] = ann.getChild(0).toXMLString().replace('&gt;', '>')
            #lists in the annotation
            #elif ann.getName()=='selenzyme' or ann.getName()=='rule_ori_reac':
            elif ann.getName()=='selenzyme':
                toRet[ann.getName()] = {}
                for y in range(ann.getNumChildren()):
                    selAnn = ann.getChild(y)
                    try:
                        toRet[ann.getName()][selAnn.getName()] = float(selAnn.getAttrValue('value'))
                    except ValueError:
                        toRet[ann.getName()][selAnn.getName()] = selAnn.getAttrValue('value')
            else:
                toRet[ann.getName()] = ann.getChild(0).toXMLString()
        #to delete empty
        return {k: v for k, v in toRet.items() if v is not None}
        #return toRet


    def readReactionSpecies(self, reaction):
        """Return the products and the species associated with a reaction

        :param reaction: Reaction object of libSBML

        :type annot: libsbml.Reaction

        :rtype: dict
        :return: Dictionary of the reaction stoichiometry
        """
        #TODO: check that reaction is either an sbml species; if not check that its a string and that
        # it exists in the rpsbml model
        toRet = {'left': {}, 'right': {}}
        #reactants
        for i in range(reaction.getNumReactants()):
            reactant_ref = reaction.getReactant(i)
            toRet['left'][reactant_ref.getSpecies()] = int(reactant_ref.getStoichiometry())
        #products
        for i in range(reaction.getNumProducts()):
            product_ref = reaction.getProduct(i)
            toRet['right'][product_ref.getSpecies()] = int(product_ref.getStoichiometry())
        return toRet


    ######################### INQUIRE ###################################


    # TODO: change the way you check for the ID since these might not be the convention
    def speciesExists(self, speciesName, compartment_id='MNXC3'):
        """Determine if the model already contains a species according to its ID

        :param reaction: Reaction object of libSBML

        :type annot: libsbml.Reaction

        :rtype: bool
        :return: True if exists and False if not
        """
        if speciesName in [i.getName() for i in self.model.getListOfSpecies()] or speciesName+'__64__'+compartment_id in [i.getId() for i in self.model.getListOfSpecies()]:
            return True
        return False


    def isSpeciesProduct(self, species_id, ignoreReactions=[]):
        """Function to determine if a species can be a product of any reaction.

        :param species_id: ID of the species to find
        :param ignoreReactions: List of all the reaction id's to ignore

        :type species_id: str
        :type ignoreReactions: list

        :rtype: bool
        :return: True if its a product of a reaction False if not
        """
        #return all the parameters values
        param_dict = {i.getId(): i.getValue() for i in self.model.parameters}
        for reaction in self.model.getListOfReactions():
            if reaction.getId() not in ignoreReactions:
                #check that the function is reversible by reversibility and FBC bounds
                if reaction.reversible:
                    reaction_fbc = reaction.getPlugin('fbc')
                    #strict left to right
                    if param_dict[reaction_fbc.getLowerFluxBound()]>=0 and param_dict[reaction_fbc.getUpperFluxBound()]>0:
                        if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                            return True
                    #can go both ways
                    elif param_dict[reaction_fbc.getLowerFluxBound()]<0 and param_dict[reaction_fbc.getUpperFluxBound()]>0:
                        if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                            return True
                        elif species_id in [i.getSpecies() for i in reaction.getListOfReactants()]:
                            return True
                    #strict right to left
                    elif param_dict[reaction_fbc.getLowerFluxBound()]<0 and param_dict[reaction_fbc.getUpperFluxBound()]<=0 and param_dict[reaction_fbc.getLowerFluxBound()]<param_dict[reaction_fbc.getUpperFluxBound()]:
                        if species_id in [i.getSpecies() for i in reaction.getListOfReactants()]:
                            return True
                    else:
                        self.logger.warning('isSpeciesProduct does not find the directionailty of the reaction for reaction: '+str(species_id))
                        return True
                else:
                    #if the reaction is not reversible then product are the only way to create it
                    if species_id in [i.getSpecies() for i in reaction.getListOfProducts()]:
                        return True
        return False


    ################### CONVERT BETWEEEN FORMATS ############################


    def outPathsDict(self, pathway_id='rp_pathway'):
        """Function to return in a dictionary in the same format as the out_paths rp2paths file dictionary object

        Example format returned: {'rule_id': 'RR-01-503dbb54cf91-49-F', 'right': {'TARGET_0000000001': 1}, 'left': {'MNXM2': 1, 'MNXM376': 1}, 'pathway_id': 1, 'step': 1, 'sub_step': 1, 'transformation_id': 'TRS_0_0_17'}. Really used to complete the monocomponent reactions

        :param pathway_id: The pathway ID (Default: rp_pathway)

        :type pathway_id: str

        :rtype: dict
        :return: Dictionary of the pathway
        """
        pathway = {}
        for member in self.getGroupsMembers(pathway_id):
            #TODO: need to find a better way
            reaction = self.model.getReaction(member)
            brsynthAnnot = self.readBRSYNTHAnnotation(reaction.getAnnotation())
            speciesReac = self.readReactionSpecies(reaction)
            self.logger.debug('brsynthAnnot:'+str(brsynthAnnot))
            step = {'reaction_id': member,
                    'reaction_rule': brsynthAnnot['smiles'],
                    'rule_score': brsynthAnnot['rule_score'],
                    'rule_id': brsynthAnnot['rule_id'],
                    'rule_ori_reac': brsynthAnnot['rule_ori_reac'],
                    'right': speciesReac['right'],
                    'left': speciesReac['left'],
                    'path_id': brsynthAnnot['path_id'],
                    'step': brsynthAnnot['step_id'],
                    'sub_step': brsynthAnnot['sub_step_id']}
            self.logger.debug('Step: '+str(step))
            pathway[brsynthAnnot['step_id']['value']] = step
        return pathway


    ############################# COMPARE MODELS ############################


    def compareBRSYNTHAnnotations(self, source_annot, target_annot):
        """Determine if two libsbml species or reactions have members in common in BRSynth annotation

        Compare two dictionnaries and if any of the values of any of the same keys are the same then the function return True, and if none are found then return False

        :param source_annot: Source object of libSBML
        :param target_annot: Target object of libSBML

        :type source_annot: libsbml.Reaction
        :type target_annot: libsbml.Reaction

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        source_dict = self.readBRSYNTHAnnotation(source_annot)
        target_dict = self.readBRSYNTHAnnotation(target_annot)
        #ignore thse when comparing reactions
        for i in ['path_id', 'step', 'sub_step', 'rule_score', 'rule_ori_reac']:
            try:
                del source_dict[i]
            except KeyError:
                pass
            try:
                del target_dict[i]
            except KeyError:
                pass
        #list the common keys between the two
        for same_key in list(set(list(source_dict.keys())).intersection(list(target_dict.keys()))):
            if source_dict[same_key]==target_dict[same_key]:
                return True
        return False


    def compareMIRIAMAnnotations(self, source_annot, target_annot):
        """Determine if two libsbml species or reactions have members in common in MIRIAM annotation

        Compare two dictionnaries and if any of the values of any of the same keys are the same then the function return True, and if none are found then return False

        :param source_annot: Source object of libSBML
        :param target_annot: Target object of libSBML

        :type source_annot: libsbml.Reaction
        :type target_annot: libsbml.Reaction

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        source_dict = self.readMIRIAMAnnotation(source_annot)
        target_dict = self.readMIRIAMAnnotation(target_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    def compareAnnotations_annot_dict(self, source_annot, target_dict):
        """Compare an annotation object and annotation dictionary

        :param source_annot: Source object of libSBML
        :param target_annot: Target dictionary

        :type target_annot: dict
        :type source_annot: libsbml.Reaction

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        source_dict = self.readMIRIAMAnnotation(source_annot)
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    def compareAnnotations_dict_dict(self, source_dict, target_dict):
        """Compare an annotation as dictionaries

        :param source_annot: Source dictionary
        :param target_annot: Target dictionary

        :type source_annot: dict
        :type target_annot: dict

        :rtype: bool
        :return: True if there is at least one similar and False if none
        """
        #list the common keys between the two
        for com_key in set(list(source_dict.keys()))-(set(list(source_dict.keys()))-set(list(target_dict.keys()))):
            #compare the keys and if same is non-empty means that there
            #are at least one instance of the key that is the same
            if bool(set(source_dict[com_key]) & set(target_dict[com_key])):
                return True
        return False


    def compareRPpathways(self, measured_sbml, pathway_id='rp_pathway'):
        """Function to compare two SBML's RP pathways

        Function that compares the annotations of reactions and if not found, the annotations of all
        species in that reaction to try to recover the correct ones. Because we are working with
        intermediate cofactors for the RP generated pathways, the annotation crossreference will
        not work. Best is to use the cross-reference to the original reaction

        :param measured_sbml: rpSBML object

        :type measured_sbml: rpSBML

        :rtype: bool, dict
        :return: True if there is at least one similar and return the dict of similarities and False if none with empty dictionary
        """
        #return all the species annotations of the RP pathways
        try:
            meas_rp_species = measured_sbml.readRPspecies()
            found_meas_rp_species = measured_sbml.readRPspecies()
            for meas_step_id in meas_rp_species:
                meas_rp_species[meas_step_id]['annotation'] = measured_sbml.model.getReaction(meas_step_id).getAnnotation()
                found_meas_rp_species[meas_step_id]['found'] = False
                for spe_name in meas_rp_species[meas_step_id]['reactants']:
                    meas_rp_species[meas_step_id]['reactants'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
                    found_meas_rp_species[meas_step_id]['reactants'][spe_name] = False
                for spe_name in meas_rp_species[meas_step_id]['products']:
                    meas_rp_species[meas_step_id]['products'][spe_name] = measured_sbml.model.getSpecies(spe_name).getAnnotation()
                    found_meas_rp_species[meas_step_id]['products'][spe_name] = False
            rp_rp_species = self.readRPspecies()
            for rp_step_id in rp_rp_species:
                rp_rp_species[rp_step_id]['annotation'] = self.model.getReaction(rp_step_id).getAnnotation()
                for spe_name in rp_rp_species[rp_step_id]['reactants']:
                    rp_rp_species[rp_step_id]['reactants'][spe_name] = self.model.getSpecies(spe_name).getAnnotation()
                for spe_name in rp_rp_species[rp_step_id]['products']:
                    rp_rp_species[rp_step_id]['products'][spe_name] = self.model.getSpecies(spe_name).getAnnotation()
        except AttributeError:
            self.logger.error('TODO: debug, for some reason some are passed as None here')
            return False, {}
        #compare the number of steps in the pathway
        if not len(meas_rp_species)==len(rp_rp_species):
            self.logger.warning('The pathways are not of the same length')
            return False, {}
        ############## compare using the reactions ###################
        for meas_step_id in measured_sbml.getGroupsMembers(pathway_id):
            for rp_step_id in rp_rp_species:
                if self.compareMIRIAMAnnotations(rp_rp_species[rp_step_id]['annotation'], meas_rp_species[meas_step_id]['annotation']):
                    found_meas_rp_species[meas_step_id]['found'] = True
                    found_meas_rp_species[meas_step_id]['rp_step_id'] = rp_step_id
                    break
        ############## compare using the species ###################
        for meas_step_id in measured_sbml.getGroupsMembers(pathway_id):
            #if not found_meas_rp_species[meas_step_id]['found']:
            for rp_step_id in rp_rp_species:
                # We test to see if the meas reaction elements all exist in rp reaction and not the opposite
                #because the measured pathways may not contain all the elements
                ########## reactants ##########
                for meas_spe_id in meas_rp_species[meas_step_id]['reactants']:
                    for rp_spe_id in rp_rp_species[rp_step_id]['reactants']:
                        if self.compareMIRIAMAnnotations(meas_rp_species[meas_step_id]['reactants'][meas_spe_id], rp_rp_species[rp_step_id]['reactants'][rp_spe_id]):
                            found_meas_rp_species[meas_step_id]['reactants'][meas_spe_id] = True
                            break
                        else:
                            if self.compareBRSYNTHAnnotations(meas_rp_species[meas_step_id]['reactants'][meas_spe_id], rp_rp_species[rp_step_id]['reactants'][rp_spe_id]):
                                found_meas_rp_species[meas_step_id]['reactants'][meas_spe_id] = True
                                break
                ########### products ###########
                for meas_spe_id in meas_rp_species[meas_step_id]['products']:
                    for rp_spe_id in rp_rp_species[rp_step_id]['products']:
                        if self.compareMIRIAMAnnotations(meas_rp_species[meas_step_id]['products'][meas_spe_id], rp_rp_species[rp_step_id]['products'][rp_spe_id]):
                            found_meas_rp_species[meas_step_id]['products'][meas_spe_id] = True
                            break
                        else:
                            if self.compareBRSYNTHAnnotations(meas_rp_species[meas_step_id]['products'][meas_spe_id], rp_rp_species[rp_step_id]['products'][rp_spe_id]):
                                found_meas_rp_species[meas_step_id]['products'][meas_spe_id] = True
                                break
                ######### test to see the difference
                pro_found = [found_meas_rp_species[meas_step_id]['products'][i] for i in found_meas_rp_species[meas_step_id]['products']]
                rea_found = [found_meas_rp_species[meas_step_id]['reactants'][i] for i in found_meas_rp_species[meas_step_id]['reactants']]
                if pro_found and rea_found:
                    if all(pro_found) and all(rea_found):
                        found_meas_rp_species[meas_step_id]['found'] = True
                        found_meas_rp_species[meas_step_id]['rp_step_id'] = rp_step_id
                        break
        ################# Now see if all steps have been found ############
        if all(found_meas_rp_species[i]['found'] for i in found_meas_rp_species):
            found_meas_rp_species['measured_model_id'] = measured_sbml.model.getId()
            found_meas_rp_species['rp_model_id'] = self.model.getId()
            return True, found_meas_rp_species
        else:
            return False, {}


    ############################# MODEL APPEND ##############################


    def setReactionConstraints(self,
                               reaction_id,
                               upper_bound,
                               lower_bound,
                               unit='mmol_per_gDW_per_hr',
                               is_constant=True):
        """Set a given reaction's upper and lower bounds

        Sets the upper and lower bounds of a reaction. Note that if the numerical values passed
        are not recognised, new parameters are created for each of them

        :param reaction_id: The id of the reaction
        :param upper_bound: Reaction upper bound
        :param lower_bound: Reaction lower bound
        :param unit: Unit to the bounds (Default: mmol_per_gDW_per_hr)
        :param is_constant: Set if the reaction is constant (Default: True)

        :type reaction_id: str
        :type upper_bound: float
        :type lower_bound: float
        :type unit: str
        :type is_constant: bool

        :rtype: tuple or bool
        :return: bool if there is an error and tuple of the lower and upper bound
        """
        reaction = self.model.getReaction(reaction_id)
        if not reaction:
            self.logger.error('Cannot find the reaction: '+str(reaction_id))
            return False
        reac_fbc = reaction.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        ########## upper bound #############
        old_upper_value = self.model.getParameter(reac_fbc.getUpperFluxBound()).value
        upper_param = self.createReturnFluxParameter(upper_bound, unit, is_constant)
        self._checklibSBML(reac_fbc.setUpperFluxBound(upper_param.getId()),
            'setting '+str(reaction_id)+' upper flux bound')
        ######### lower bound #############
        old_lower_value = self.model.getParameter(reac_fbc.getLowerFluxBound()).value
        lower_param = self.createReturnFluxParameter(lower_bound, unit, is_constant)
        self._checklibSBML(reac_fbc.setLowerFluxBound(lower_param.getId()),
            'setting '+str(reaction_id)+' lower flux bound')
        return old_upper_value, old_lower_value


    ##### ADD SOURCE FROM ORPHAN #####
    #if the heterologous pathway from the self.model contains a sink molecule that is not included in the 
    # original model (we call orhpan species) then add another reaction that creates it
    #TODO: that transports the reactions that creates the species in the
    # extracellular matrix and another reaction that transports it from the extracellular matrix to the cytoplasm
    #TODO: does not work
    def fillOrphan(self,
            rpsbml=None,
            pathway_id='rp_pathway',
            compartment_id='MNXC3',
            upper_flux_bound=999999,
            lower_flux_bound=10):
        """Fill the orgpan

        WARNING: in progress

        :rtype: tuple or bool
        :return: bool if there is an error and tuple of the lower and upper bound
        """
        if rpsbml==None:
            model = self.model
        else:
            model = rpsbml.model
        self.logger.debug('Adding the orphan species to the GEM model')
        #only for rp species
        groups = model.getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        reaction_id = sorted([(int(''.join(x for x in i.id_ref if x.isdigit())), i.id_ref) for i in rp_pathway.getListOfMembers()], key=lambda tup: tup[0], reverse=True)[0][1]
        #for reaction_id in [i.getId() for i in self.model.getListOfReactions()]:
        for species_id in set([i.getSpecies() for i in model.getReaction(reaction_id).getListOfReactants()]+[i.getSpecies() for i in model.getReaction(reaction_id).getListOfProducts()]):
            if rpsbml==None:
                isSpePro = self.isSpeciesProduct(species_id, [reaction_id])
            else:
                isSpePro = rpsbml.isSpeciesProduct(species_id, [reaction_id])
            if not isSpePro:
                #create the step
                createStep = {'rule_id': None,
                              'left': {species_id.split('__')[0]: 1},
                              'right': {},
                              'step': None,
                              'sub_step': None,
                              'path_id': None,
                              'transformation_id': None,
                              'rule_score': None,
                              'rule_ori_reac': None}
                #create the model in the
                if rpsbml==None:
                    self.createReaction('create_'+species_id,
                                        upper_flux_bound,
                                        lower_flux_bound,
                                        createStep,
                                        compartment_id)
                else:
                    rpsbml.createReaction('create_'+species_id,
                                        upper_flux_bound,
                                        lower_flux_bound,
                                        createStep,
                                        compartment_id)


    ############################# MODEL CREATION FUNCTIONS ##################


    def createModel(self, name, model_id, meta_id=None):
        """Create libSBML model instance

        Function that creates a new libSBML model instance and initiates it with the appropriate packages. Creates a cytosol compartment

        :param name: The name of the of the model
        :param model_id: The id of the model
        :param meta_id: Meta ID of the model (Default: None)

        :type name: str
        :type model_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        ## sbmldoc
        self.sbmlns = libsbml.SBMLNamespaces(3,1)
        self._checklibSBML(self.sbmlns, 'generating model namespace')
        self._checklibSBML(self.sbmlns.addPkgNamespace('groups',1), 'Add groups package')
        self._checklibSBML(self.sbmlns.addPkgNamespace('fbc',2), 'Add FBC package')
        #sbmlns = libsbml.SBMLNamespaces(3,1,'groups',1)
        self.document = libsbml.SBMLDocument(self.sbmlns)
        self._checklibSBML(self.document, 'generating model doc')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('fbc', False), 'enabling FBC package')
        #!!!! must be set to false for no apparent reason
        self._checklibSBML(self.document.setPackageRequired('groups', False), 'enabling groups package')
        ## sbml model
        self.model = self.document.createModel()
        self._checklibSBML(self.model, 'generating the model')
        self._checklibSBML(self.model.setId(model_id), 'setting the model ID')
        model_fbc = self.model.getPlugin('fbc')
        model_fbc.setStrict(True)
        if meta_id==None:
            meta_id = self._genMetaID(model_id)
        self._checklibSBML(self.model.setMetaId(meta_id), 'setting model meta_id')
        self._checklibSBML(self.model.setName(name), 'setting model name')
        self._checklibSBML(self.model.setTimeUnits('second'), 'setting model time unit')
        self._checklibSBML(self.model.setExtentUnits('mole'), 'setting model compartment unit')
        self._checklibSBML(self.model.setSubstanceUnits('mole'), 'setting model substance unit')


    #TODO: set the compName as None by default. To do that you need to regenerate the compXref to
    #TODO: consider seperating it in another function if another compartment is to be created
    #TODO: use MNX ids as keys instead of the string names
    def createCompartment(self, size, compId, compName, compXref, meta_id=None):
        """Create libSBML compartment

        :param size: Size of the compartment
        :param compId: Compartment id
        :param compName: Compartment Name
        :param compXref: Cross reference dictionary of the compartment
        :param meta_id: Meta id (Default: None)
        :type size: float

        :type compId: str
        :type compName: str
        :type compXref: dict
        :type meta_id: str

        :rtype: None
        :return: None
        """
        comp = self.model.createCompartment()
        self._checklibSBML(comp, 'create compartment')
        self._checklibSBML(comp.setId(compId), 'set compartment id')
        if compName:
            self._checklibSBML(comp.setName(compName), 'set the name for the cytoplam')
        self._checklibSBML(comp.setConstant(True), 'set compartment "constant"')
        self._checklibSBML(comp.setSize(size), 'set compartment "size"')
        self._checklibSBML(comp.setSBOTerm(290), 'set SBO term for the cytoplasm compartment')
        if meta_id==None:
            meta_id = self._genMetaID(compId)
        self._checklibSBML(comp.setMetaId(meta_id), 'set the meta_id for the compartment')
        ############################ MIRIAM ############################
        comp.setAnnotation(libsbml.XMLNode.convertStringToXMLNode(self._defaultMIRIAMAnnot(meta_id)))
        self.addUpdateMIRIAM(comp, 'compartment', compXref, meta_id)


    def createUnitDefinition(self, unit_id, meta_id=None):
        """Create libSBML unit definition

        Function that creates a unit definition (composed of one or more units)

        :param unit_id: Unit id definition
        :param meta_id: Meta id (Default: None)

        :type unit_id: str
        :type meta_id: str

        :rtype: libsbml.UnitDefinition
        :return: Unit definition object created
        """
        unitDef = self.model.createUnitDefinition()
        self._checklibSBML(unitDef, 'creating unit definition')
        self._checklibSBML(unitDef.setId(unit_id), 'setting id')
        if meta_id==None:
            meta_id = self._genMetaID(unit_id)
        self._checklibSBML(unitDef.setMetaId(meta_id), 'setting meta_id')
        #self.unitDefinitions.append(unit_id)
        return unitDef


    def createUnit(self, unitDef, libsbmlunit, exponent, scale, multiplier):
        """Set or update the parameters of a libSBML unit definition

        :param unitDef: libSBML Unit
        :param libsbmlunit: String unit
        :param exponent: Exponent unit
        :param sale: Scale of the unit
        :param multiplier: Multiplier of the unit

        :type unitDef: libsbml.Unit
        :type libsbmlunit: str
        :type exponent: int
        :type sale: int
        :type multiplier: int

        :rtype: None
        :return: None
        """
        unit = unitDef.createUnit()
        self._checklibSBML(unit, 'creating unit')
        self._checklibSBML(unit.setKind(libsbmlunit), 'setting the kind of unit')
        self._checklibSBML(unit.setExponent(exponent), 'setting the exponenent of the unit')
        self._checklibSBML(unit.setScale(scale), 'setting the scale of the unit')
        self._checklibSBML(unit.setMultiplier(multiplier), 'setting the multiplier of the unit')


    def createReturnFluxParameter(self,
                                  value,
                                  unit='mmol_per_gDW_per_hr',
                                  is_constant=True,
                                  parameter_id=None,
                                  meta_id=None):
        """Create libSBML flux parameters

        Parameters are used for the bounds for FBA analysis. Unit parameter must be an instance of unitDefinition.
        If the parameter id exists, then the function returns the libsbml.Parameter object

        :param value: Value set for the parameter
        :param unit: The unit id of the parameter (Default: mmol_per_gDW_per_hr)
        :param is_constant: Define if the parameter is constant (Default: True)
        :param parameter_id: Overwrite the default naming convention (Default: None)
        :param meta_id: Meta id (Default: None)

        :type value: float
        :type unit: str
        :type is_constant: bool
        :type parameter_id: str
        :type meta_id: str

        :rtype: libsbml.Parameter
        :return: The newly created libsbml.Parameter
        """
        if parameter_id:
            param_id = parameter_id
        else:
            if value>=0:
                param_id = 'B_'+str(round(abs(value), 4)).replace('.', '_')
            else:
                param_id = 'B__'+str(round(abs(value), 4)).replace('.', '_')
        if param_id in [i.getId() for i in self.model.getListOfParameters()]:
            return self.model.getParameter(param_id)
        else:
            newParam = self.model.createParameter()
            self._checklibSBML(newParam, 'Creating a new parameter object')
            self._checklibSBML(newParam.setConstant(is_constant), 'setting as constant')
            self._checklibSBML(newParam.setId(param_id), 'setting ID')
            self._checklibSBML(newParam.setValue(value), 'setting value')
            self._checklibSBML(newParam.setUnits(unit), 'setting units')
            self._checklibSBML(newParam.setSBOTerm(625), 'setting SBO term')
            if meta_id==None:
                meta_id = self._genMetaID(parameter_id)
            self._checklibSBML(newParam.setMetaId(meta_id), 'setting meta ID')
            #self.parameters.append(parameter_id)
            return newParam


    #TODO as of now not generic, works when creating a new SBML file, but no checks if modifying existing SBML file
    def createReaction(self,
                       reac_id,
                       fluxUpperBound,
                       fluxLowerBound,
                       step,
                       compartment_id,
                       reaction_smiles=None,
                       reacXref={},
                       pathway_id=None,
                       meta_id=None):
        """Create libSBML reaction

        Create a reaction that is added to the self.model in the input compartment id. fluxBounds is a list of libSBML.UnitDefinition, length of exactly 2 with the first position that is the upper bound and the second is the lower bound. reactants_dict and reactants_dict are dictionnaries that hold the following parameters: name, compartment, stoichiometry

        :param name: Name of the reaction
        :param fluxUpperBound: The reaction fbc upper bound
        :param fluxLowerBound: The reaction fbc lower bound
        :param step: The id's of the reactant and products of the reactions. Example: {'left': [], 'right': []}
        :param compartment_id: The id of the compartment to add the reaction
        :param reaction_smiles: The reaction rule to add to the BRSynth annotation of the reaction (Default: None)
        :param reacXref: The dict containing the MIRIAM annotation (Default: {})
        :param pathway_id: The Groups id of the reaction to which the reacion id will be added (Default: None)
        :param meta_id: Meta id (Default: None)

        :type name: str
        :type fluxUpperBound: float
        :type fluxLowerBound: float
        :type step: dict
        :type compartment_id: str
        :type reaction_smiles: str
        :type reacXref: dict
        :type pathway_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        reac = self.model.createReaction()
        self._checklibSBML(reac, 'create reaction')
        ################ FBC ####################
        reac_fbc = reac.getPlugin('fbc')
        self._checklibSBML(reac_fbc, 'extending reaction for FBC')
        #bounds
        upper_bound = self.createReturnFluxParameter(fluxUpperBound)
        lower_bound = self.createReturnFluxParameter(fluxLowerBound)
        self._checklibSBML(reac_fbc.setUpperFluxBound(upper_bound.getId()), 'setting '+str(reac_id)+' upper flux bound')
        self._checklibSBML(reac_fbc.setLowerFluxBound(lower_bound.getId()), 'setting '+str(reac_id)+' lower flux bound')
        #########################################
        #reactions
        self._checklibSBML(reac.setId(reac_id), 'set reaction id') #same convention as cobrapy
        self._checklibSBML(reac.setSBOTerm(176), 'setting the system biology ontology (SBO)') #set as process
        #TODO: consider having the two parameters as input to the function
        self._checklibSBML(reac.setReversible(True), 'set reaction reversibility flag')
        self._checklibSBML(reac.setFast(False), 'set reaction "fast" attribute')
        if meta_id==None:
            meta_id = self._genMetaID(reac_id)
        self._checklibSBML(reac.setMetaId(meta_id), 'setting species meta_id')
        #TODO: check that the species exist
        #reactants_dict
        for reactant in step['left']:
            spe = reac.createReactant()
            self._checklibSBML(spe, 'create reactant')
            #use the same writing convention as CobraPy
            self._checklibSBML(spe.setSpecies(str(reactant)+'__64__'+str(compartment_id)), 'assign reactant species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(spe.setConstant(True), 'set "constant" on species '+str(reactant))
            self._checklibSBML(spe.setStoichiometry(float(step['left'][reactant])),
                'set stoichiometry ('+str(float(step['left'][reactant]))+')')
        #TODO: check that the species exist
        #products_dict
        for product in step['right']:
            pro = reac.createProduct()
            self._checklibSBML(pro, 'create product')
            self._checklibSBML(pro.setSpecies(str(product)+'__64__'+str(compartment_id)), 'assign product species')
            #TODO: check to see the consequences of heterologous parameters not being constant
            self._checklibSBML(pro.setConstant(True), 'set "constant" on species '+str(product))
            self._checklibSBML(pro.setStoichiometry(float(step['right'][product])),
                'set the stoichiometry ('+str(float(step['right'][product]))+')')
        ############################ MIRIAM ############################
        self._checklibSBML(reac.setAnnotation(self._defaultBothAnnot(meta_id)), 'creating annotation')
        self.addUpdateMIRIAM(reac, 'reaction', reacXref, meta_id)
        ###### BRSYNTH additional information ########
        if reaction_smiles:
            self.addUpdateBRSynth(reac, 'smiles', reaction_smiles, None, True, False, False, meta_id)
        if step['rule_id']:
            self.addUpdateBRSynth(reac, 'rule_id', step['rule_id'], None, True, False, False, meta_id)
        #TODO: need to change the name and content (to dict) upstream
        if step['rule_ori_reac']:
            self.addUpdateBRSynth(reac, 'rule_ori_reac', step['rule_ori_reac'], None, True, False, False, meta_id)
        if step['rule_score']:
            self.addUpdateBRSynth(reac, 'rule_score', step['rule_score'], None, False, False, False, meta_id)
        if step['path_id']:
            self.addUpdateBRSynth(reac, 'path_id', step['path_id'], None, False, False, False, meta_id)
        if step['step']:
            self.addUpdateBRSynth(reac, 'step_id', step['step'], None, False, False, False, meta_id)
        if step['sub_step']:
            self.addUpdateBRSynth(reac, 'sub_step_id', step['sub_step'], None, False, False, False, meta_id)
        #### GROUPS #####
        if not pathway_id==None:
            groups_plugin = self.model.getPlugin('groups')
            hetero_group = groups_plugin.getGroup(pathway_id)
            if not hetero_group:
                self.logger.warning('The pathway_id '+str(pathway_id)+' does not exist in the model')
            else:
                newM = hetero_group.createMember()
                self._checklibSBML(newM, 'Creating a new groups member')
                self._checklibSBML(newM.setIdRef(reac_id), 'Setting name to the groups member')


    def createSpecies(self,
                      species_id,
                      compartment_id,
                      species_name=None,
                      chemXref={},
                      inchi=None,
                      inchikey=None,
                      smiles=None,
                      species_group_id=None,
                      in_sink_group_id=None,
                      meta_id=None):
                      #TODO: add these at some point -- not very important
                      #charge=0,
                      #chemForm=''):
        """Create libSBML species

        Create a species that is added to self.model

        :param species_id: The id of the created species
        :param compartment_id: The id of the compartment to add the reaction
        :param species_name: Overwrite the default name of the created species (Default: None)
        :param chemXref: The dict containing the MIRIAM annotation (Default: {})
        :param inchi: The InChI string to be added to BRSynth annotation (Default: None)
        :param inchikey: The InChIkey string to be added to BRSynth annotation (Default: None)
        :param smiles: The SMLIES string to be added to BRSynth annotation (Default: None)
        :param species_group_id: The Groups id to add the species (Default: None)
        :param in_sink_group_id: The Groups id sink species to add the species (Default: None)
        :param meta_id: Meta id (Default: None)

        :type species_id: str
        :type compartment_id: str
        :type species_name: str
        :type chemXref: dict
        :type inchi: str
        :type inchikey: str
        :type smiles: str
        :type species_group_id: str
        :type in_sink_group_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        spe = self.model.createSpecies()
        self._checklibSBML(spe, 'create species')
        ##### FBC #####
        spe_fbc = spe.getPlugin('fbc')
        self._checklibSBML(spe_fbc, 'creating this species as an instance of FBC')
        #spe_fbc.setCharge(charge) #### These are not required for FBA
        #spe_fbc.setChemicalFormula(chemForm) #### These are not required for FBA
        #if compartment_id:
        self._checklibSBML(spe.setCompartment(compartment_id), 'set species spe compartment')
        #else:
        #    #removing this could lead to errors with xref
        #    self._checklibSBML(spe.setCompartment(self.compartment_id), 'set species spe compartment')
        #ID same structure as cobrapy
        #TODO: determine if this is always the case or it will change
        self._checklibSBML(spe.setHasOnlySubstanceUnits(False), 'set substance units')
        self._checklibSBML(spe.setBoundaryCondition(False), 'set boundary conditions')
        self._checklibSBML(spe.setConstant(False), 'set constant')
        #useless for FBA (usefull for ODE) but makes Copasi stop complaining
        self._checklibSBML(spe.setInitialConcentration(1.0), 'set an initial concentration')
        #same writting convention as COBRApy
        self._checklibSBML(spe.setId(str(species_id)+'__64__'+str(compartment_id)), 'set species id')
        self.logger.debug('Setting species id as: '+str(species_id)+'__64__'+str(compartment_id))
        if meta_id==None:
            meta_id = self._genMetaID(species_id)
        self._checklibSBML(spe.setMetaId(meta_id), 'setting reaction meta_id')
        if species_name==None:
            self._checklibSBML(spe.setName(species_id), 'setting name for the metabolite '+str(species_id))
        else:
            self._checklibSBML(spe.setName(species_name), 'setting name for the metabolite '+str(species_name))
        #this is setting MNX id as the name
        #this is setting the name as the input name
        #self._checklibSBML(spe.setAnnotation(self._defaultBRSynthAnnot(meta_id)), 'creating annotation')
        self._checklibSBML(spe.setAnnotation(self._defaultBothAnnot(meta_id)), 'creating annotation')
        ###### annotation ###
        self.addUpdateMIRIAM(spe, 'species', chemXref, meta_id)
        ###### BRSYNTH additional information ########
        if smiles:
            self.addUpdateBRSynth(spe, 'smiles', smiles, None, True, False, False, meta_id)
            #                   sbase_obj, annot_header, value, units=None, isAlone=False, isList=False, isSort=True, meta_id=None)
        if inchi:
            self.addUpdateBRSynth(spe, 'inchi', inchi, None, True, False, False, meta_id)
        if inchikey:
            self.addUpdateBRSynth(spe, 'inchikey', inchikey, None, True, False, False, meta_id)
            self.addUpdateMIRIAM(spe, 'species', {'inchikey': [inchikey]})
        #### GROUPS #####
        #TODO: check that it actually exists
        self.logger.debug('species_group_id: '+str(species_group_id))
        if not species_group_id==None:
            groups_plugin = self.model.getPlugin('groups')
            hetero_group = groups_plugin.getGroup(species_group_id)
            if not hetero_group:
                self.logger.warning('The species_group_id '+str(species_group_id)+' does not exist in the model')
                #TODO: consider creating it if
            else:
                newM = hetero_group.createMember()
                self._checklibSBML(newM, 'Creating a new groups member')
                self._checklibSBML(newM.setIdRef(str(species_id)+'__64__'+str(compartment_id)), 'Setting name to the groups member')
        #TODO: check that it actually exists
        #add the species to the sink species
        self.logger.debug('in_sink_group_id: '+str(in_sink_group_id))
        if not in_sink_group_id==None:
            groups_plugin = self.model.getPlugin('groups')
            sink_group = groups_plugin.getGroup(in_sink_group_id)
            if not sink_group:
                self.logger.warning('The species_group_id '+str(in_sink_group_id)+' does not exist in the model')
                #TODO: consider creating it if
            else:
                newM = sink_group.createMember()
                self._checklibSBML(newM, 'Creating a new groups member')
                self._checklibSBML(newM.setIdRef(str(species_id)+'__64__'+str(compartment_id)), 'Setting name to the groups member')


    #TODO: change the name of this function to createGroup
    def createPathway(self, pathway_id, meta_id=None):
        """Create libSBML pathway

        Create a pathway that is added to self.model

        :param pathway_id: The Groups id of the pathway id
        :param meta_id: Meta id (Default: None)

        :type pathway_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        groups_plugin = self.model.getPlugin('groups')
        new_group = groups_plugin.createGroup()
        self.logger.debug('setting new group id: '+str(pathway_id))
        new_group.setId(pathway_id)
        if meta_id==None:
            meta_id = self._genMetaID(pathway_id)
        new_group.setMetaId(meta_id)
        new_group.setKind(libsbml.GROUP_KIND_COLLECTION)
        new_group.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    def createGene(self, reac, step_id, meta_id=None):
        """Create libSBML gene

        Create a gene that is associated with a reaction

        :param reac: The id of the reaction that is associated with the gene
        :param step_id: The id of the reaction to name the gene
        :param meta_id: Meta id (Default: None)

        :type reac: str
        :type step_id: str
        :type meta_id: str

        :rtype: None
        :return: None
        """
        #TODO: pass this function to Pablo for him to fill with parameters that are appropriate for his needs
        geneName = 'RP'+str(step_id)+'_gene'
        fbc_plugin = self.model.getPlugin('fbc')
        #fbc_plugin = reac.getPlugin("fbc")
        gp = fbc_plugin.createGeneProduct()
        gp.setId(geneName)
        if meta_id==None:
            meta_id = self._genMetaID(str(geneName))
        gp.setMetaId(meta_id)
        gp.setLabel('gene_'+str(step_id))
        gp.setAssociatedSpecies('RP'+str(step_id))
        ##### NOTE: The parameters here require the input from Pablo to determine what he needs
        #gp.setAnnotation(self._defaultBothAnnot(meta_id))


    def createFluxObj(self, fluxobj_id, reaction_name, coefficient, is_max=True, meta_id=None):
        """Create libSBML flux objective

        WARNING DEPRECATED -- use the createMultiFluxObj() with lists of size one to define an objective function
        with a single reaction
        Using the FBC package one can add the FBA flux objective directly to the model. This function sets a particular reaction as objective with maximization or minimization objectives

        :param fluxobj_id: The id of the flux objective
        :param reaction_name: The id of the reaction that is associated with the reaction
        :param coefficient: The coefficient of the flux objective
        :param is_max: Define if the objective is coefficient (Default: True)
        :param meta_id: Meta id (Default: None)

        :type fluxobj_id: str
        :type reaction_name: str
        :type coefficient: int
        :type is_max: bool
        :type meta_id: str

        :rtype: None
        :return: None
        """
        fbc_plugin = self.model.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        #TODO: need to define inpiut metaID
        target_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        target_obj.setId(fluxobj_id)
        if is_max:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxobj_id) # this ensures that we are using this objective when multiple
        target_flux_obj = target_obj.createFluxObjective()
        target_flux_obj.setReaction(reaction_name)
        target_flux_obj.setCoefficient(coefficient)
        if meta_id==None:
            meta_id = self._genMetaID(str(fluxobj_id))
        target_flux_obj.setMetaId(meta_id)
        target_flux_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    def createMultiFluxObj(self, fluxobj_id, reaction_names, coefficients, is_max=True, meta_id=None):
        """Create libSBML flux objective

        Using the FBC package one can add the FBA flux objective directly to the model. Can add multiple reactions. This function sets a particular reaction as objective with maximization or minimization objectives

        :param fluxobj_id: The id of the flux objective
        :param reaction_names: The list of string id's of the reaction that is associated with the reaction
        :param coefficients: The list of int defining the coefficients of the flux objective
        :param is_max: Define if the objective is coefficient (Default: True)
        :param meta_id: Meta id (Default: None)

        :type fluxobj_id: str
        :type reaction_names: list
        :type coefficients: list
        :type is_max: bool
        :type meta_id: str

        :rtype: None
        :return: None
        """
        if not len(reaction_names)==len(coefficients):
            self.logger.error('The size of reaction_names is not the same as coefficients')
            return False
        fbc_plugin = self.model.getPlugin('fbc')
        target_obj = fbc_plugin.createObjective()
        target_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))
        target_obj.setId(fluxobj_id)
        if is_max:
            target_obj.setType('maximize')
        else:
            target_obj.setType('minimize')
        fbc_plugin.setActiveObjectiveId(fluxobj_id) # this ensures that we are using this objective when multiple
        for reac, coef in zip(reaction_names, coefficients):
            target_flux_obj = target_obj.createFluxObjective()
            target_flux_obj.setReaction(reac)
            target_flux_obj.setCoefficient(coef)
            if meta_id==None:
                meta_id = self._genMetaID(str(fluxobj_id))
            target_flux_obj.setMetaId(meta_id)
            target_flux_obj.setAnnotation(self._defaultBRSynthAnnot(meta_id))


    ##############################################################################################
    ############################### Generic Model ################################################
    ##############################################################################################


    def genericModel(self,
                     model_name,
                     model_id,
                     compXref,
                     compartment_id,
                     upper_flux_bound=999999,
                     lower_flux_bound=0):
        """Generate a generic model

        Since we will be using the same type of parameters for the RetroPath model, this function
        generates a libSBML model with parameters that will be mostly used

        :param model_name: The given name of the model
        :param model_id: The id of the model
        :param compXref: The model MIRIAM annotation
        :param compartment_id: The id of the model compartment
        :param upper_flux_bound: The upper flux bounds unit definitions default when adding new reaction (Default: 999999.0)
        :param lower_flux_bound: The lower flux bounds unit definitions default when adding new reaction (Defaul: 0.0)

        :type model_name: str
        :type model_id: str
        :type compXref: dict
        :type compartment_id: str
        :type upper_flux_bound: float
        :type lower_flux_bound: float

        :rtype: None
        :return: None
        """
        self.createModel(model_name, model_id)
        # mmol_per_gDW_per_hr -- flux
        unitDef = self.createUnitDefinition('mmol_per_gDW_per_hr')
        self.createUnit(unitDef, libsbml.UNIT_KIND_MOLE, 1, -3, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_GRAM, 1, 0, 1)
        self.createUnit(unitDef, libsbml.UNIT_KIND_SECOND, 1, 0, 3600)
        # kj_per_mol -- thermodynamics
        gibbsDef = self.createUnitDefinition('kj_per_mol')
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_JOULE, 1, 3, 1)
        self.createUnit(gibbsDef, libsbml.UNIT_KIND_MOLE, -1, 1, 1)
        ### set the bounds
        upBound = self.createReturnFluxParameter(upper_flux_bound)
        lowBound = self.createReturnFluxParameter(lower_flux_bound)
        #compartment
        #TODO: try to recover the name from the Xref
        try:
            name = compXref['name'][0]
        except KeyError:
            name = compartment_id+'_name'
        self.createCompartment(1, compartment_id, name, compXref)

