import csv
import requests
import itertools
import re
import time
import os
import json
from io import StringIO
import copy
import logging

from .rpSBML import rpSBML
from .rpCache import rpCache

__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = ["Joan Herisson"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Melchior du Lac"
__status__ = "Development"


logging.basicConfig(
    level=logging.DEBUG,
    #level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


class Species:
    """Class holding all the possible gathered information that we have for a chemical species
    """
    def __init__(self, cid, inchi=None, inchikey=None, smiles=None, name=None, xref=None):
        """The constructor
        """
        self.inchi = inchi
        self.inchikey = inchikey
        self.smiles = smiles
        self.xref = xref
        self.name = name
        self.cid = cid
    def __eq__(self, cid):
        return self.cid==cid
    def __hash__(self):
        return hash(self.cid)


class rpReader(rpCache):
    """Class that has all the functions that parse different files to convert them to rpSBML objects and files

    Supports: - RetroPath2.0 input
              - TSV
              - String input

    .. document private functions
    .. automethod:: _pubChemLimit
    .. automethod:: _pubchemStrctSearch
    """
    def __init__(self, rpcache=None):
        """Constructor for the class
        """
        #we do this for speed - although given inheritance you would think that the data is loaded twice, the MEM addresses are passed here
        super().__init__()
        if rpcache:
            super().__init__()
            self.deprecatedCID_cid = rpcache.deprecatedCID_cid
            self.deprecatedRID_rid = rpcache.deprecatedRID_rid
            self.cid_strc = rpcache.cid_strc
            self.cid_xref = rpcache.cid_xref
            self.comp_xref = rpcache.comp_xref
            self.xref_comp = rpcache.xref_comp
            self.cid_name = rpcache.cid_name
            self.chebi_cid = rpcache.chebi_cid
            self.inchikey_cid = rpcache.inchikey_cid
            self.rr_reactions = rpcache.rr_reactions
            self.rr_full_reactions = rpcache.rr_full_reactions
        self.rpcache = rpcache
        self.logger = logging.getLogger(__name__)
        self.logger.info('Starting instance of rpReader')
        #species
        self.species_obj_dict = []


    #######################################################################
    ######################### Rp2paths and RetroPath2 #####################
    #######################################################################


    def rp2ToSBML(self,
                  rp2_pathways,
                  rp2paths_compounds,
                  rp2paths_pathways,
                  out_dir,
                  upper_flux_bound=999999.0,
                  lower_flux_bound=0.0,
                  max_subpaths_filter=100,
                  pathway_id='rp_pathway',
                  compartment_id='MNXC3',
                  species_group_id='central_species',
                  sink_species_group_id='rp_sink_species',
                  pubchem_search=False):
        """Function to group all the functions for parsing RP2 output to SBML files

        Takes RP2paths's compounds.txt and out_paths.csv and RetroPaths's *_scope.csv files and generates SBML

        :param rp2_pathways: The RetroPath2.0 results scope file
        :param rp2paths_pathways: The rp2paths result pathway (out_paths) file
        :param rp2paths_compounds: The rp2paths result compounds file
        :param tmpOutputFolder: A folder to output the results (Default: None)
        :param upper_flux_bound: The default upper flux bound (Default: 999999)
        :param lower_flux_bound: The default lower flux bound (Default: 0)
        :param max_subpaths_filter: The maximal number of rules associated with each step (Default: 10)
        :param pathway_id: The Groups heterologous pathway id (Default: rp_pathway)
        :param compartment_id: The compartment SBML id (Default: MNXC3)
        :param species_group_id: The Groups id of the central species (Default: central_species)
        :param sink_species_group_id: The Groups id of the rp_sink_species (Default: rp_sink_species)
        :param pubchem_search: Use the pubchem database to search for missing cross reference (Default: False)

        :type rp2_pathways: str
        :type rp2paths_pathways: str
        :type rp2paths_compounds: str
        :type tmpOutputFolder: str
        :type upper_flux_bound: int
        :type lower_flux_bound: int
        :type max_subpaths_filter: int
        :type pathway_id: str
        :type compartment_id: str
        :type species_group_id: str
        :type sink_species_group_id: str
        :type pubchem_search: bool

        :rtype: dict
        :return: Dictionnary of pathways results
        """
        if max_subpaths_filter<0:
            raise ValueError('Max number of subpaths cannot be less than 0: '+str(max_subpaths_filter))
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        rp_strc = self.compounds(rp2paths_compounds)
        rp_transformation, sink_molecules = self.transformations(rp2_pathways)
        return self.rp2pathsToSBML(rp_strc,
                                   rp_transformation,
                                   sink_molecules,
                                   rp2paths_pathways,
                                   out_dir,
                                   upper_flux_bound,
                                   lower_flux_bound,
                                   max_subpaths_filter,
                                   pathway_id,
                                   compartment_id,
                                   species_group_id,
                                   sink_species_group_id,
                                   pubchem_search)


    def compounds(self, path):
        """Function to parse the compounds.txt file

        Extract the smile and the structure of each compounds of RP2Path output. Method to parse all the RP output compounds.

        :param path: Path to the compounds file

        :type path: str

        :rtype: dict
        :return: Dictionnary of compounds results
        """
        rp_strc = {}
        try:
            if isinstance(path, bytes):
                reader = csv.reader(StringIO(path.decode('utf-8')), delimiter='\t')
            else:
                reader = csv.reader(open(path, 'r', encoding='utf-8'), delimiter='\t')
            next(reader)
            for row in reader:
                rp_strc[row[0]] = {'smiles': row[1]}  #, 'structure':row[1].replace('[','').replace(']','')
                strc_info = None
                try:
                    if not strc_info:
                        strc_info = self.queryCIDstrc(row[0])
                    rp_strc[row[0]]['inchi'] = strc_info['inchi']
                except KeyError:
                    #try to generate them yourself by converting them directly
                    try:
                        resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchi'})
                        rp_strc[row[0]]['inchi'] = resConv['inchi']
                    except NotImplementedError as e:
                        logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                try:
                    if not strc_info:
                        strc_info = self.queryCIDstrc(row[0])
                    rp_strc[row[0]]['inchikey'] = strc_info['inchikey']
                    #try to generate them yourself by converting them directly
                    #TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
                except KeyError:
                    try:
                        resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})
                        rp_strc[row[0]]['inchikey'] = resConv['inchikey']
                    except NotImplementedError as e:
                        logger.warning('Could not convert the following SMILES to InChI key: '+str(row[1]))
        except (TypeError, FileNotFoundError) as e:
            logger.error('Could not read the compounds file ('+str(path)+')')
            raise RuntimeError
        return rp_strc


    def transformations(self, path):
        """Function to parse the scope.csv file

        Extract the reaction rules from the retroPath2.0 output using the scope.csv file

        :param path: Path to the compounds file

        :type path: str

        :rtype: tuple
        :return: The RetroPath transformation and the list of sink molecules
        """
        rp_transformation = {}
        sink_molecules = []
        #### we might pass binary in the REST version
        reader = None
        if isinstance(path, bytes):
            reader = csv.reader(StringIO(path.decode('utf-8')), delimiter=',')
        else:
            try:
                reader = csv.reader(open(path, 'r'), delimiter=',')
            except FileNotFoundError:
                logger.error('Could not read the compounds file: '+str(path))
                return {}
        next(reader)
        for row in reader:
            if not row[1] in rp_transformation:
                rp_transformation[row[1]] = {}
                rp_transformation[row[1]]['rule'] = row[2]
                rp_transformation[row[1]]['ec'] = [i.replace(' ', '') for i in row[11][1:-1].split(',') if not i.replace(' ', '')=='NOEC']
            if row[7]=='1':
                for i in row[8].replace(']', '').replace('[', '').replace(' ', '').split(','):
                    sink_molecules.append(i)
        # logger.info(rp_transformation)
        # logger.info(sink_molecules)
        return rp_transformation, list(set(sink_molecules))


    def readRp2Paths(self, rp2paths_pathways):
        """Function that reads the pathway output of rp2paths

        :param rp2paths_pathways: The path to the rp2paths pathway output

        :type rp2paths_pathways: str

        :rtype: dict
        :return: Return dictionary of rp2paths
        """
        #### we might pass binary in the REST version
        if isinstance(rp2paths_pathways, bytes):
            reader = csv.reader(StringIO(rp2paths_pathways.decode('utf-8')))
        else:
            reader = csv.reader(open(rp2paths_pathways, 'r'))
        next(reader)
        current_path_id = 0
        path_step = 1
        rp_paths = {}
        for row in reader:
            try:
                #Remove all illegal characters in SBML ids
                row[3] = row[3].replace("'", "").replace('-', '_').replace('+', '')
                if not int(row[0])==current_path_id:
                    path_step = 1
                else:
                    path_step += 1
                #important to leave them in order
                current_path_id = int(row[0])
            except ValueError:
                logger.error('Cannot convert path_id to int ('+str(row[0])+')')
                #return {}
                return False
            #################################
            ruleIds = row[2].split(',')
            if ruleIds==None:
                logger.warning('The rulesIds is None')
                #pass # or continue
                continue
            ###WARNING: This is the part where we select some rules over others
            # we do it by sorting the list according to their score and taking the topx
            tmp_rr_reactions = {}
            for r_id in ruleIds:
                for rea_id in self.rr_reactions[r_id]:
                    tmp_rr_reactions[str(r_id)+'__'+str(rea_id)] = self.rr_reactions[r_id][rea_id]
            # if len(ruleIds)>int(max_subpaths_filter):
            #     logger.warning('There are too many rules, limiting the number to random top '+str(max_subpaths_filter))
            #     try:
            #         ruleIds = [y for y,_ in sorted([(i, tmp_rr_reactions[i]['rule_score']) for i in tmp_rr_reactions])][:int(max_subpaths_filter)]
            #     except KeyError:
            #         logger.warning('Could not select topX due inconsistencies between rules ids and rr_reactions... selecting random instead')
            #         ruleIds = random.sample(tmp_rr_reactions, int(max_subpaths_filter))
            # else:
            ruleIds = tmp_rr_reactions
            sub_path_step = 1
            for singleRule in ruleIds:
                tmpReac = {'rule_id': singleRule.split('__')[0],
                           'rule_ori_reac': singleRule.split('__')[1],
                           'rule_score': self.rr_reactions[singleRule.split('__')[0]][singleRule.split('__')[1]]['rule_score'],
                           'right': {},
                           'left': {},
                           'path_id': int(row[0]),
                           'step': path_step,
                           'transformation_id': row[1][:-2]}
                ############ LEFT ##############
                for l in row[3].split(':'):
                    tmp_l = l.split('.')
                    #tmpReac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                    cid = '' #TODO: change this
                    cid = self._checkCIDdeprecated(tmp_l[1])
                    try:
                        tmpReac['left'][cid] = int(tmp_l[0])
                    except ValueError:
                        logger.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                        #return {}
                        return False
                ############## RIGHT ###########
                for r in row[4].split(':'):
                    tmp_r = r.split('.')
                    #tmpReac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                    cid = '' #TODO change this
                    cid = self._checkCIDdeprecated(tmp_r[1])
                    try:
                        tmpReac['right'][cid] = int(tmp_r[0])
                    except ValueError:
                        logger.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                        return False
                #################################
                if not int(row[0]) in rp_paths:
                    rp_paths[int(row[0])] = {}
                if not int(path_step) in rp_paths[int(row[0])]:
                    rp_paths[int(row[0])][int(path_step)] = {}
                rp_paths[int(row[0])][int(path_step)][int(sub_path_step)] = tmpReac
                sub_path_step += 1
        return rp_paths


    def _gatherSpeciesInfo(self, species, pubchem_search=False):
        """Given a Species object, try to fill all the missing information

        :param species: Species object
        :param pubchem_search: Query the pubchem webb service

        :type species: Species
        :type pubchem_search: bool

        :rtype: bool
        :return: Success or failure of the function
        """
        #check that the cid is valid
        species.cid = self._checkCIDdeprecated(species.cid)
        #if there are no cid and but the inchikey, try to recover 
        if not species.cid in self.cid_xref or (not species.cid and species.inchikey):
            try:
                species.cid = self._checkCIDdeprecated([i for i in self.inchikey_cid[species.inchikey] if i[:3]=='MNX'][0])
            except KeyError:
                self.logger.error('Cannot find the cid using the inchikey: '+str(species.inchikey))
                return False
        if not species.cid:
            self.logger.error('There must be a cid of the Species object')
            return False
        #common name of the species
        if not species.name:
            try:
                #species.name = self.cid_strc[species.cid]['name'].replace("'", "")
                species.name = self.queryCIDname(species.cid)
            except KeyError:
                species.name = None
        ###Try to fill from cache
        if not species.xref:
            try:
                species.xref = self.cid_xref[species.cid]
            except KeyError:
                species.xref = {}
        otype = set()
        in_strct = None
        itype = None
        strc_info = None
        if not species.inchi:
            try:
                if not strc_info:
                    strc_info = self.queryCIDstrc(species.cid)
                species.inchi = strc_info['inchi']
                if not in_strct:
                    in_strct = species.inchi
                    itype = 'inchi'
            except KeyError:
                otype.add('inchi')
        if not species.smiles:
            try:
                if not strc_info:
                    strc_info = self.queryCIDstrc(species.cid)
                species.smiles = strc_info['smiles']
                if not in_strct:
                    in_strct = species.smiles
                    itype = 'smiles'
            except KeyError:
                otype.add('smiles')
        if not species.inchikey:
            try:
                if not strc_info:
                    strc_info = self.queryCIDstrc(species.cid)
                species.inchikey = strc_info['inchikey']
                if not in_strct:
                    in_strct = species.inchikey
                    itype = 'inchikey'
            except KeyError:
                otype.add('inchikey')
        ###Try with pubchem
        if ((in_strct and itype) or not species.xref) and pubchem_search:
            try:
                pubres = self._pubchemStrctSearch(in_strct, itype)
                if not chem_name:
                    chem_name = pubres['name']
                if 'chebi' in pubres['xref']:
                    try:
                        species.xref = self.cid_xref[self.chebi_cid[pubres['xref']['chebi'][0]]]
                    except KeyError:
                        pass
                if not species.xref:
                    species.xref = pubres['xref']
                if not species.inchikey:
                    species.inchikey = pubres['inchikey']
                    if 'inchikey' in otype:
                        otype.remove('inchikey')
                if not species.smiles:
                    species.smiles = pubres['smiles']
                    if 'smiles' in otype:
                        otype.remove('smiles')
                if not species.inchi:
                    species.inchi = pubres['inchi']
                    if 'inchi' in otype:
                        otype.remove('inchi')
            except KeyError:
                self.logger.warning('Bad results from pubchem results for query: '+str(species.inchi))
        ###Try to convert from known information
        if otype and itype and in_strct:
            strc_conv = self._convert_depiction(in_strc, itype, otype)
            for s in strc_conv:
                if s=='inchi' and not species.inchi:
                    species.inchi = strc_conv[s]
                elif s=='inchikey' and not species.inchikey:
                    species.inchikey = strc_conv[s]
                elif s=='smiles' and not species.smiles:
                    species.smiles = strc_conv[s]
        return True


    def rp2pathsToSBML(self,
                       rp_strc,
                       rp_transformation,
                       sink_molecules,
                       rp2paths_pathways,
                       out_folder,
                       upper_flux_bound=999999.0,
                       lower_flux_bound=0.0,
                       max_subpaths_filter=10,
                       pathway_id='rp_pathway',
                       compartment_id='MNXC3',
                       species_group_id='central_species',
                       sink_species_group_id='rp_sink_species',
                       pubchem_search=False):
        """Function to convert the already parsed rp2paths file to SBML files

        Reading the RP2path output and extract all the information for each pathway RP2path Metabolic pathways from out_paths.csv create all the different values for heterologous paths from the RP2path out_paths.csv file. Note that path_step are in reverse order here

        :param rp_strc: The structures dictionary from RetroPath2.0
        :param rp_transformation: The transformation dictionary from RetroPath2.0
        :param sink_molecules: The dictionary of sink molecules from RetroPath2.0
        :param rp2paths_pathways: The rp2paths result file
        :param out_folder: A folder to output the results (Default: None)
        :param upper_flux_bound: The default upper flux bound (Default: 999999)
        :param lower_flux_bound: The default lower flux bound (Default: 0)
        :param max_subpaths_filter: The maximal number of rules associated with each step (Default: 10)
        :param pathway_id: The Groups heterologous pathway id (Default: rp_pathway)
        :param compartment_id: The compartment SBML id (Default: MNXC3)
        :param species_group_id: The Groups id of the central species (Default: central_species)
        :param sink_species_group_id: The Groups id of the rp_sink_species (Default: rp_sink_species)
        :param pubchem_search: Use the pubchem database to search for missing cross reference (Default: False)

        :type rp_strc: str
        :type rp_transformation: str
        :type sink_molecules: str
        :type rp2paths_pathways: str
        :type out_folder: str
        :type upper_flux_bound: int
        :type lower_flux_bound: int
        :type max_subpaths_filter: int
        :type pathway_id: str
        :type compartment_id: str
        :type species_group_id: str
        :type sink_species_group_id: str
        :type pubchem_search: bool

        :rtype: dict
        :return: The dictionary of the rp2paths results or False boolean if fails
        """
        # TODO: make sure that you account for the fact that each reaction may have multiple associated reactions
        rp_paths = self.readRp2Paths(rp2paths_pathways)
        sink_species = []
        # for each line or rp2paths_pathways:
        #     generate comb
        #     for each combinant:
        #         rank
        #         process
        #         add cofactors
        #         dedup
        #### pathToSBML ####
        try:
            compartment_id = self.xref_comp[compartment_id]
        except KeyError:
            logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
            return False
        for pathNum in rp_paths:
            # first level is the list of lists of sub_steps
            # second is itertools all possible combinations using product
            alt_path_num = 1
            # sopX subpaths of the current rp2path pathway
            list_rpsbml = []
            for comb_path in list(itertools.product(*[[(i,y) for y in rp_paths[pathNum][i]] for i in rp_paths[pathNum]])):
                steps = []
                for i, y in comb_path:
                    steps.append(rp_paths[pathNum][i][y])
                path_id = steps[0]['path_id']
                rpsbml = rpSBML(model_name='rp_'+str(path_id)+'_'+str(alt_path_num), rpcache=self.rpcache)
                #1) Create a generic Model, ie the structure and unit definitions that we will use the most
                ##### TODO: give the user more control over a generic model creation:
                # -> special attention to the compartment
                rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(alt_path_num),
                                    'RP_model_'+str(path_id)+'_'+str(alt_path_num),
                                    compartment_id,
                                    upper_flux_bound=upper_flux_bound,
                                    lower_flux_bound=lower_flux_bound)
                #2) Create the pathway (groups)
                rpsbml.createGroup(pathway_id)
                rpsbml.createGroup(species_group_id)
                rpsbml.createGroup(sink_species_group_id)
                #3) Find all unique species and add them to the model
                all_cid = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
                for cid in all_cid:
                    cid = self._checkCIDdeprecated(cid)
                    species = Species(cid=cid)
                    # inchi
                    try:
                        species.inchi = rp_strc[cid]['inchi']
                    except KeyError:
                        pass
                    #inchikey
                    try:
                        species.inchikey = rp_strc[cid]['inchikey']
                    except KeyError:
                        pass
                    #smiles
                    try:
                        species.smiles = rp_strc[cid]['smiles']
                    except KeyError:
                        pass
                    if species in self.species_obj_dict:
                        species = self.species_obj_dict[self.species_obj_dict.index(species)]
                    else:
                        self._gatherSpeciesInfo(species, pubchem_search)
                    #pass the information to create the species
                    #rpsbml = self._addSpecies(rpsbml, cid, sink_molecules, compartment_id, chem_name, spe, species_group_id, sink_species_group_id)
                    #create the species according to the type of group they belong to
                    if not species in self.species_obj_dict:
                        self.species_obj_dict.append(species)
                    if cid in sink_molecules:
                        rpsbml.createSpecies(cid,
                                             compartment_id,
                                             species.name,
                                             species.xref,
                                             species.inchi,
                                             species.inchikey,
                                             species.smiles,
                                             species_group_id,
                                             sink_species_group_id)
                    else:
                        rpsbml.createSpecies(cid,
                                             compartment_id,
                                             species.name,
                                             species.xref,
                                             species.inchi,
                                             species.inchikey,
                                             species.smiles,
                                             species_group_id)
                #4) Add the complete reactions and their annotations
                for step in steps:
                    # add the substep to the model
                    step['sub_step'] = alt_path_num
                    rpsbml.createReaction('RP'+str(step['step']), # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['step']),
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          step,
                                          compartment_id,
                                          rp_transformation[step['transformation_id']]['rule'],
                                          {'ec': rp_transformation[step['transformation_id']]['ec']},
                                          pathway_id)
                #5) Adding the consumption of the target
                target_step = {'rule_id': None,
                               'left': {[i for i in all_cid if i[:6]=='TARGET'][0]: 1}, #warning this is dangerous
                               'right': {},
                               'step': None,
                               'sub_step': None,
                               'path_id': None,
                               'transformation_id': None,
                               'rule_score': None,
                               'rule_ori_reac': None}
                rpsbml.createReaction('RP1_sink',
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      target_step,
                                      compartment_id)
                #6) Adding the cofactors
                self._addCofactors(rpsbml, compartment_id, pathway_id, pubchem_search)
                #7) Insert the new rpsbml object in list if it does not exist
                #NOTE: sorting or bisect is a O^N problem, so no difference
                if not rpsbml in list_rpsbml:
                    list_rpsbml.append(rpsbml)
                    list_rpsbml.sort(reverse=True)
                    # 8) Keep only topX
                    list_rpsbml = list_rpsbml[:max_subpaths_filter]
                alt_path_num += 1
            # Write results to files
            for rpsbml in list_rpsbml:
                self.logger.debug('out_folder: '+str(out_folder))
                rpsbml.writeSBML(out_folder)
        return True


    #############################################################################################
    ############################### COFACTORS ###################################################
    #############################################################################################


    def completeReac(self, step, rr_reac, full_reac, mono_side, rr_string, pathway_cmp):
        """Given a dictionary describing a monocomponent reaction, add the cofactors by comparing it with the original reaction

        :param step: Dictionnary describing the reaction
        :param rr_reac: Dictionnary describing the monocomponent reaction from RetroRules
        :param full_reac: The original full reaction description
        :param mono_side: Is monocomponent side of the reaction
        :param rr_string: The reaction rule
        :param pathway_cmp: Dictionnary used to retreive the public ID of the intermediate compounds. Resets for each individual pathway

        :type step: dict
        :type rr_reac: dict
        :type full_reac: dict
        :type mono_side: dict
        :type rr_string: str
        :type pathway_cmp: dict

        :rtype: tuple
        :return: The tuple with the status of the function and the complete reaction string
        """
        if mono_side:
            ## add the unknown species to pathway_cmp for the next steps
            rr_mono_cmp = list(rr_reac.keys())
            step_mono_cmp = list(step.keys())
            if (len(rr_mono_cmp)==1 and len(step_mono_cmp)==1):
                #this is purposely overwitten since the main cmp between reactions can change
                pathway_cmp[step_mono_cmp[0]] = rr_mono_cmp[0]
            else:
                logger.warning('There should be only one compound on the left for monocomponent reaction: rr_mono_cmp: '+str(rr_mono_cmp)+' step_mono_cmp: '+str(step_mono_cmp))
                return False
        ## add the side species
        rr_string += self._addCofactorSpecies(step, full_reac, rr_reac)
        ## Update the stochio
        return True, self._updateStoichio(step, full_reac, rr_string, pathway_cmp)


    def _addCofactorSpecies(self, step, full_reac, rr_reac):
        """Add the new cofactor species and update the reaction rule

        :param step: The dictionary describing the reaction
        :param full_reac: The original full reaction
        :param rr_reac: The Reaction rule reaction describing

        :type step: dict
        :type full_reac: dict
        :type rr_reac: dict

        :rtype: str
        :return: The updated reaction rule string with the cofactors
        """
        rr_string = ''
        for toAdd in full_reac.keys()-rr_reac.keys():
            step.update({toAdd: full_reac[toAdd]})
            ### update the reaction rule string
            try:
                smi = self.cid_strc[toAdd]['smiles']
                if not smi==None:
                    for sto_add in range(int(full_reac[toAdd])):
                        rr_string += '.'+str(smi)
            except KeyError:
                self.logger.warning('Cannot find smiles structure for '+str(toAdd))
        return rr_string


    def _updateStoichio(self, step, full_reac, rr_string, pathway_cmp):
        """Update the stoichiometry of the reactions

        :param step: The dictionary describing the reaction
        :param full_reac: The original full reaction
        :param rr_string: The reaction rule
        :param pathway_cmp: The compounds in the pathway

        :type step: dict
        :type full_reac: dict
        :type rr_string: str
        :type pathway_cmp: dict

        :rtype: str
        :return: The updated reaction rule
        """
        for step_spe in step:
            if step_spe in full_reac:
                if not step[step_spe]==full_reac[step_spe]:
                    stochio_diff = full_reac[step_spe]-step[step_spe]
                    step[step_spe] = full_reac[step_spe]
                    if stochio_diff<0:
                        logger.warning('full_reac stochio should never be smaller than step')
                        continue
                    for i in range(stochio_diff):
                        ### update the reaction rule string
                        try:
                            smi = self.cid_strc[step_spe]['smiles']
                            if not smi==None:
                                rr_string += '.'+str(smi)
                        except KeyError:
                            #@Mel toAdd -> step_spe
                            logger.warning('Cannot find smiles structure for '+str(step_spe))
            elif step_spe in pathway_cmp:
                if pathway_cmp[step_spe] in full_reac:
                    if not step[step_spe]==full_reac[pathway_cmp[step_spe]]:
                        step[step_spe] = full_reac[pathway_cmp[step_spe]]
            #Its fine if the stochio is not updated, better than ignoring a whole pathway
                #else:
                #    logger.warning('Cannot find '+str(step_spe)+' in full reaction')
                #    return False
            #else:
            #    logger.warning('Cannot find '+str(step_spe)+' in pathway_cmp')
            #    return False
            return rr_string


    def _addCofactorsStep(self, step, pathway_cmp):
        """Add the cofactors to monocomponent reactions

        :param step: Step in a pathway
        :param pathway_cmp: Intermediate compounds with their public ID's

        :type step: dict
        :type pathway_cmp: dict

        :rtype: bool
        :return: Success or failure of the function
        """
        reac_smiles_left = step['reaction_rule'].split('>>')[0]
        reac_smiles_right = step['reaction_rule'].split('>>')[1]
        if self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==-1:
            try:
                isSuccess, reac_smiles_left = self.completeReac(step['right'],
                        self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                        self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'])]['right'],
                        True,
                        reac_smiles_left,
                        pathway_cmp)
                if not isSuccess:
                    logger.warning('Could not recognise reaction rule for step (1): '+str(step))
                    return False
            except KeyError:
                logger.warning('Could not find the full reaction for reaction (1): '+str(step))
                return False
            try:
                isSuccess, reac_smiles_right = self.completeReac(step['left'],
                        self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                        self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'])]['left'],
                        False,
                        reac_smiles_right,
                        pathway_cmp)
                if not isSuccess:
                    logger.warning('Could not recognise reaction rule for step (2): '+str(step))
                    return False
            except KeyError:
                logger.warning('Could not find the full reaction for reaction (2): '+str(step))
                return False
        elif self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==1:
            try:
                isSuccess, reac_smiles_left = self.completeReac(step['right'],
                        self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                        self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'])]['left'],
                        True,
                        reac_smiles_left,
                        pathway_cmp)
                if not isSuccess:
                    logger.error('Could not recognise reaction rule for step (3): '+str(step))
                    return False
            except KeyError:
                logger.warning('Could not find the full reaction for reaction (3): '+str(step))
                return False
            try:
                isSuccess, reac_smiles_right = self.completeReac(step['left'],
                        self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                        self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'])]['right'],
                        False,
                        reac_smiles_right,
                        pathway_cmp)
                if not isSuccess:
                    logger.error('Could not recognise reaction rule for step (4): '+str(step))
                    return False
            except KeyError:
                logger.warning('Could not find the full reaction for reaction (4): '+str(step))
                return False
        else:
            logger.error('Relative direction can only be 1 or -1: '+str(self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']))
            return False
        step['reaction_rule'] = reac_smiles_left+'>>'+reac_smiles_right
        return True


    def _addCofactors(self, rpsbml, compartment_id='MNXC3', pathway_id='rp_pathway', pubchem_search=False):
        """Function to reconstruct the heterologous pathway

        Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors

        :param rpsbml: The rpSBML object with a single model
        :param compartment_id: The id of the SBML compartment of interest
        :param pathway_id: The Groups id of the heterologous pathway
        :param pubchem_search: Query the pubchem database

        :type rpsbml: rpSBML
        :type compartment_id: str
        :type pathway_id: str
        :type pubchem_search: bool

        :rtype: bool
        :return: Success or failure of the function
        """
        #This keeps the IDs conversions to the pathway
        pathway_cmp = {}
        spe_conv = {}
        rpsbml_dict = rpsbml.asdict(pathway_id)
        rp_path = rpsbml.outPathsDict(pathway_id)
        ori_rp_path = copy.deepcopy(rp_path)
        #We reverse the loop to ID the intermediate CMP to their original ones
        for step_num in sorted(list(rp_path), reverse=True):
        #for step_num in sorted(list(rp_path)):
            if self._addCofactorsStep(rp_path[step_num], pathway_cmp):
                ###add the new cofactors to the SBML
                #remove the original species from the monocomponent reaction
                reactants = set(set(rp_path[step_num]['left'].keys())-set(ori_rp_path[step_num]['left'].keys()))
                products = set(set(rp_path[step_num]['right'].keys())-set(ori_rp_path[step_num]['right'].keys()))
                for in_cid in reactants|products:
                    self.logger.debug('in_cid: '+str(in_cid))
                    cid = self._checkCIDdeprecated(in_cid)
                    self.logger.debug('in_cid: '+str(in_cid))
                    self.logger.debug('cid: '+str(cid))
                    #check to make sure that they do not yet exist and if not create a new one
                    #TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                    if not rpsbml.speciesExists(cid, compartment_id):
                        #check if the cid exists 
                        try:
                            species = Species(cid=cid, inchikey=rpsbml_dict['species'][in_cid]['brsynth']['inchikey'])
                        except KeyError:
                            species = Species(cid=cid)
                        if species in self.species_obj_dict:
                            species = self.species_obj_dict[self.species_obj_dict.index(species)]
                        else:
                            self._gatherSpeciesInfo(species, pubchem_search)
                        #### Create the species in the SBML file ######
                        rpsbml.createSpecies(species.cid,
                                             compartment_id,
                                             species.name,
                                             species.xref,
                                             species.inchi,
                                             species.inchikey,
                                             species.smiles)
                        if not species in self.species_obj_dict:
                            self.species_obj_dict.append(species)
                #add the new species to the RP reactions
                reac = rpsbml.model.getReaction(rp_path[step_num]['reaction_id'])
                if not reac:
                    self.logger.error('Cannot find the following reaction: '+str(rp_path[step_num]['reaction_id']))
                    return False
                pre_reactants = [i.species for i in reac.getListOfReactants()]
                pre_products = [i.species for i in reac.getListOfProducts()]
                for pro in products:
                    if self._checkCIDdeprecated(pro) in spe_conv:
                        toadd = spe_conv[self._checkCIDdeprecated(pro)]
                    else:
                        toadd = str(self._checkCIDdeprecated(pro))+'__64__'+str(compartment_id)
                    #prod.setSpecies(str(self._checkCIDdeprecated(pro))+'__64__'+str(compartment_id))
                    if toadd in pre_products:
                        continue
                    prod = reac.createProduct()
                    prod.setSpecies(toadd)
                    prod.setConstant(True)
                    prod.setStoichiometry(rp_path[step_num]['right'][pro])
                for sub in reactants:
                    if self._checkCIDdeprecated(sub) in spe_conv:
                        toadd = spe_conv[self._checkCIDdeprecated(sub)]
                    else:
                        toadd = str(self._checkCIDdeprecated(sub))+'__64__'+str(compartment_id)
                    #prod.setSpecies(str(self._checkCIDdeprecated(sub))+'__64__'+str(compartment_id))
                    if toadd in pre_reactants:
                        continue
                    subs = reac.createReactant()
                    subs.setSpecies(toadd)
                    subs.setConstant(True)
                    subs.setStoichiometry(rp_path[step_num]['left'][sub])
                #replace the reaction rule with new one
                rpsbml.addUpdateBRSynth(reac, 'smiles', rp_path[step_num]['reaction_rule'], None, True)
            else:
                #if the cofactors cannot be found delete it from the list
                logger.warning('Cannot find cofactors... skipping')
                return False
        return True


    #############################################################################################
    ######################################## TSV ################################################
    #############################################################################################

    #TODO: need to update to Species

    ## Function to parse the TSV of measured heterologous pathways to SBML
    #
    # TODO: update this to the new compartements and others
    # Given the TSV of measured pathways, parse them to a dictionary, readable to next be parsed
    # to SBML
    #
    # @param self object pointer
    # @param inFile The input JSON file
    # @param mnxHeader Reorganise the results around the target MNX products
    # @return Dictionnary of SBML
    def _parseTSV(self, inFile, remove_inchi_4p=False, mnxHeader=False):
        data = {}
        try:
            for row in csv.DictReader(open(inFile), delimiter='\t'):
                ######## path_id ######
                try:
                    pathID = int(row['pathway_ID'])
                except ValueError:
                    logger.error('Cannot convert pathway ID: '+str(row['pathway_ID']))
                    continue
                if not pathID in data:
                    data[pathID] = {}
                    data[pathID]['isValid'] = True
                    data[pathID]['steps'] = {}
                ####### target #########
                if not 'target' in data[pathID]:
                    data[pathID]['target'] = {}
                    data[pathID]['target']['name'] = row['target_name']
                    if remove_inchi_4p:
                        data[pathID]['target']['inchi'] = '/'.join([row['target_structure'].split('/')[i] for i in range(len(row['target_structure'].split('/'))) if i<4])
                    else:
                        data[pathID]['target']['inchi'] = row['target_structure']
                ####### step #########
                try:
                    stepID = int(row['step'])
                except ValueError:
                    logger.error('Cannot convert step ID: '+str(row['step']))
                    data[pathID]['isValid'] = False
                    continue
                if stepID==0:
                    continue
                elif stepID==1:
                    data[pathID]['organism'] = row['organism'].replace(' ', '')
                    data[pathID]['reference'] = row['reference'].replace(' ', '')
                data[pathID]['steps'][stepID] = {}
                ##### substrates #########
                data[pathID]['steps'][stepID]['substrates'] = []
                lenDBref = len(row['substrate_dbref'].split(';'))
                for i in row['substrate_dbref'].split(';'):
                    if i=='':
                        lenDBref -= 1
                lenStrc = len(row['substrate_structure'].split('_'))
                for i in row['substrate_structure'].split('_'):
                    if i=='':
                        lenStrc -= 1
                lenSub = len(row['substrate_name'].split(';'))
                for i in row['substrate_name'].split(';'):
                    if i=='':
                        lenSub -= 1
                if lenSub==lenStrc==lenSub:
                    for name, inchi, dbrefs in zip(row['substrate_name'].split(';'),
                            row['substrate_structure'].split('_'),
                            row['substrate_dbref'].split(';')):
                        tmp = {}
                        if remove_inchi_4p:
                            tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
                        else:
                            tmp['inchi'] = inchi.replace(' ', '')
                        tmp['name'] = name
                        tmp['dbref'] = {}
                        for dbref in dbrefs.split('|'):
                            if len(dbref.split(':'))==2:
                                db_name = dbref.split(':')[0].replace(' ', '').lower()
                                db_cid = dbref.split(':')[1].replace(' ', '')
                                if not db_name in tmp['dbref']:
                                    tmp['dbref'][db_name] = []
                                tmp['dbref'][db_name].append(db_cid)
                            else:
                                logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                                data[pathID]['isValid'] = False
                        data[pathID]['steps'][stepID]['substrates'].append(tmp)
                else:
                    logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['substrate_name'])+' <--> '+str(row['substrate_structure'])+' <--> '+str(row['substrate_dbref']))
                    data[pathID]['isValid'] = False
                    continue
                ##### products #########
                data[pathID]['steps'][stepID]['products'] = []
                lenDBref = len(row['product_dbref'].split(';'))
                for i in row['product_dbref'].split(';'):
                    if i=='':
                        lenDBref -= 1
                lenStrc = len(row['product_structure'].split('_'))
                for i in row['product_structure'].split('_'):
                    if i=='':
                        lenStrc -= 1
                lenSub = len(row['product_name'].split(';'))
                for i in row['product_name'].split(';'):
                    if i=='':
                        lenSub -= 1
                if lenSub==lenStrc==lenDBref:
                    for name, inchi, dbrefs in zip(row['product_name'].split(';'),
                            row['product_structure'].split('_'),
                            row['product_dbref'].split(';')):
                        tmp = {}
                        if remove_inchi_4p:
                            tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
                        else:
                            tmp['inchi'] = inchi.replace(' ', '')
                        tmp['name'] = name
                        tmp['dbref'] = {}
                        for dbref in dbrefs.split('|'):
                            if len(dbref.split(':'))==2:
                                db_name = dbref.split(':')[0].replace(' ', '').lower()
                                db_cid = dbref.split(':')[1].replace(' ', '')
                                if not db_name in tmp['dbref']:
                                    tmp['dbref'][db_name] = []
                                tmp['dbref'][db_name].append(db_cid)
                            else:
                                data[pathID]['isValid'] = False
                                logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                        data[pathID]['steps'][stepID]['products'].append(tmp)
                else:
                    logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['product_name'])+' <--> '+str(row['product_structure'])+' <--> '+str(row['product_dbref']))
                    data[pathID]['isValid'] = False
                if not row['uniprot']=='':
                    data[pathID]['steps'][stepID]['uniprot'] = row['uniprot'].replace(' ', '').split(';')
                if not row['EC_number']=='':
                    data[pathID]['steps'][stepID]['ec_numbers'] = [i.replace(' ', '') for i in row['EC_number'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_id'] = [i.replace(' ', '') for i in row['enzyme_identifier'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_name'] = row['enzyme_name'].split(';')
        except FileNotFoundError:
            logger.error('Cannot open the file: '+str(inFile))
        #now loop through all of them and remove the invalid paths
        toRet = copy.deepcopy(data)
        for path_id in data.keys():
            if toRet[path_id]['isValid']==False:
                del toRet[path_id]
            else:
                del toRet[path_id]['isValid']
        #reorganise the results around the target products mnx
        if not mnxHeader:
            return toRet
        else:
            toRetTwo = {}
            for path_id in toRet:
                try:
                    final_pro_mnx = toRet[path_id]['steps'][max(toRet[path_id]['steps'])]['products'][0]['dbref']['mnx'][0]
                except KeyError:
                    logger.error('The species '+str(toRet[path_id]['steps'][max(toRet[path_id]['steps'])]['products'][0]['name'])+' does not contain a mnx database reference... skipping whole pathway number '+str(path_id))
                    #continue
                if not final_pro_mnx in toRetTwo:
                    toRetTwo[final_pro_mnx] = {}
                toRetTwo[final_pro_mnx][path_id] = toRet[path_id]
            return toRetTwo


    ## Parse the validation TSV to SBML
    #
    # Parse the TSV file to SBML format and adds them to the sbml_paths
    #
    # @param self Object pointer
    # @param inFile Input file
    # @param compartment_id compartment of the
    # TODO: update this with the new SBML groups
    def TSVtoSBML(self,
                  inFile,
                  output_folder=None,
                  upper_flux_bound=99999,
                  lower_flux_bound=0.0,
                  compartment_id='MNXC3',
                  pathway_id='rp_pathway',
                  species_group_id='central_species',
                  header_name=''):
        data = self._parseTSV(inFile)
        sbml_paths = {}
        if header_name=='':
            header_name = inFile.split('/')[-1].replace('.tsv', '').replace('.csv', '')
        # TODO: need to exit at this loop
        for path_id in data:
            try:
                mnxc = self.xref_comp[compartment_id]
            except KeyError:
                logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
                return False
            rpsbml = rpSBML.rpSBML(model_name=header_name+'_'+str(path_id), rpcache=self.rpcache)
            # 1) create a generic Model, ie the structure and unit definitions that we will use the most
            ##### TODO: give the user more control over a generic model creation:
            # -> special attention to the compartment
            rpsbml.genericModel(header_name+'_Path'+str(path_id),
                                header_name+'_Path'+str(path_id),
                                self.comp_xref[mnxc],
                                compartment_id,
                                upper_flux_bound,
                                lower_flux_bound)
            # 2) create the pathway (groups)
            rpsbml.createGroup(pathway_id)
            rpsbml.createGroup(species_group_id)
            # 3) find all the unique species and add them to the model
            allChem = []
            for stepNum in data[path_id]['steps']:
                # because of the nature of the input we need to remove duplicates
                for i in data[path_id]['steps'][stepNum]['substrates']+data[path_id]['steps'][stepNum]['products']:
                    if not i in allChem:
                        allChem.append(i)
            # add them to the SBML
            for chem in allChem:
                # PROBLEM: as it stands one expects the cid to be MNX
                if 'mnx' in chem['dbref']:
                    # must list the different models
                    cid = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                else:
                    # TODO: add the species with other types of xref in annotation
                    logger.warning('Some species are not referenced by a MNX id and will be ignored')
                    # try CHEBI
                    try:
                        cid = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                        cid = 'CHEBI_'+str(cid)
                    except KeyError:
                        # TODO: need to find a better way
                        logger.warning('Cannot determine MNX or CHEBI entry, using random')
                        tmpDB_name = list(chem['dbref'].keys())[0]
                        cid = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                        cid = str(tmpDB_name)+'_'+str(cid)
                    # break
                # try to conver the inchi into the other structures
                smiles = None
                inchikey = None
                try:
                    resConv = self._convert_depiction(idepic=chem['inchi'], itype='inchi', otype={'smiles','inchikey'})
                    smiles = resConv['smiles']
                    inchikey = resConv['inchikey']
                except NotImplementedError as e:
                    logger.warning('Could not convert the following InChI: '+str(chem['inchi']))
                # create a new species
                # here we want to gather the info from rpReader's rp_strc and cid_strc
                try:
                    chem_name = self.cid_strc[cid]['name']
                except KeyError:
                    chem_name = cid
                # compile as much info as you can
                # xref
                try:
                    # TODO: add the xref from the document
                    spe_xref = self.cid_xref[cid]
                except KeyError:
                    #spe_xref = {}
                    spe_xref = chem['dbref']
                # inchi
                try:
                    spe_inchi = self.cid_strc[cid]['inchi']
                except KeyError:
                    spe_inchi = chem['inchi']
                # inchikey
                try:
                    spe_inchikey = self.cid_strc[cid]['inchikey']
                except KeyError:
                    spe_inchikey =  resConv['inchikey']
                # smiles
                try:
                    spe_smiles = self.cid_strc[cid]['smiles']
                except KeyError:
                    spe_smiles = resConv['smiles']
                # pass the information to create the species
                rpsbml.createSpecies(cid,
                                     compartment_id,
                                     chem_name,
                                     spe_xref,
                                     spe_inchi,
                                     spe_inchikey,
                                     spe_smiles,
                                     species_group_id)
            # 4) add the complete reactions and their annotations
            # create a new group for the measured pathway
            # need to convert the validation to step for reactions
            for stepNum in data[path_id]['steps']:
                toSend = {'left': {}, 'right': {}, 'rule_id': None, 'rule_ori_reac': None, 'rule_score': None, 'path_id': path_id, 'step': stepNum, 'sub_step': None}
                for chem in data[path_id]['steps'][stepNum]['substrates']:
                    if 'mnx' in chem['dbref']:
                        cid = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        # try CHEBI
                    else:
                        logger.warning('Not all the species to have a MNX ID')
                        # break
                        try:
                            cid = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            cid = 'CHEBI_'+str(cid)
                        except KeyError:
                            # TODO: need to find a better way
                            logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmpDB_name = list(chem['dbref'].keys())[0]
                            cid = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            cid = str(tmpDB_name)+'_'+str(cid)
                    toSend['left'][cid] = 1
                for chem in data[path_id]['steps'][stepNum]['products']:
                    if 'mnx' in chem['dbref']:
                        cid = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        # try CHEBI
                    else:
                        logger.warning('Need all the species to have a MNX ID')
                        try:
                            cid = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            cid = 'CHEBI_'+str(cid)
                        except KeyError:
                            # TODO: need to find a better way
                            logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmpDB_name = list(chem['dbref'].keys())[0]
                            cid = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            cid = str(tmpDB_name)+'_'+str(cid)
                    toSend['right'][cid] = 1
                        # break
                # if all are full add it
                reac_xref = {}
                if 'ec_numbers' in data[path_id]['steps'][stepNum]:
                    reac_xref['ec'] = data[path_id]['steps'][stepNum]['ec_numbers']
                if 'uniprot' in data[path_id]['steps'][stepNum]:
                    reac_xref['uniprot'] = data[path_id]['steps'][stepNum]['uniprot']
                logger.debug('#########################################')
                logger.debug(toSend)
                logger.debug('#########################################')
                rpsbml.createReaction(header_name+'_Step'+str(stepNum),
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      toSend,
                                      compartment_id,
                                      None,
                                      reac_xref,
                                      pathway_id)
                if stepNum==1:
                    # adding the consumption of the target
                    target_step = {'rule_id': None,
                                   'left': {},
                                   'right': {},
                                   'step': None,
                                   'sub_step': None,
                                   'path_id': None,
                                   'transformation_id': None,
                                   'rule_score': None,
                                   'rule_ori_reac': None}
                    for chem in data[path_id]['steps'][stepNum]['products']:
                        try:
                            # smallest MNX
                            cid = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        except KeyError:
                            # try CHEBI
                            try:
                                cid = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                                cid = 'CHEBI_'+str(cid)
                            except KeyError:
                                logger.warning('Cannot determine MNX or CHEBI entry, using random')
                                tmpDB_name = list(chem['dbref'].keys())[0]
                                cid = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                                cid = str(tmpDB_name)+'_'+str(cid)
                        target_step['left'][cid] = 1
                    rpsbml.createReaction(header_name+'_Step1_sink',
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          target_step,
                                          compartment_id)
                    rpsbml.createFluxObj('rpFBA_obj', header_name+'_Step1_sink', 1, True)
            if output_folder:
                rpsbml.writeSBML(output_folder)
            else:
                sbml_paths[header_name+'_Path'+str(path_id)] = rpsbml
        if output_folder:
            return {}
        else:
            return sbml_paths
