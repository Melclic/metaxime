import csv
import glob
import requests
import itertools
import tempfile
import tarfile
import re
import time
import os
import json
from io import StringIO
import copy
import logging
import time

from .rpSBML import rpSBML
from .rpCache import rpCache

__author__ = "Melchior du Lac"
__copyright__ = "Copyright 2020"
__credits__ = ["Joan Herisson"]
__license__ = "GPLv3"
__version__ = "0.0.1"
__maintainer__ = "Melchior du Lac"
__status__ = "Development"


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
            self.deprecatedCID_cid = rpcache.deprecatedCID_cid
            self.deprecatedRID_rid = rpcache.deprecatedRID_rid
            self.deprecatedComp_comp = rpcache.deprecatedComp_comp
            self.cid_strc = rpcache.cid_strc
            self.cid_xref = rpcache.cid_xref
            self.comp_xref = rpcache.comp_xref
            self.cid_name = rpcache.cid_name
            self.chebi_cid = rpcache.chebi_cid
            self.inchikey_cid = rpcache.inchikey_cid
            self.rr_reactions = rpcache.rr_reactions
            self.rr_full_reactions = rpcache.rr_full_reactions
        self.rpcache = rpcache
        #self.logger = logging.getLogger(__name__)
        self.logger = logging.getLogger(os.path.basename(__file__))
        self.logger.info('Starting instance of rpReader')
        #species
        self.species_obj_dict = []


    #######################################################################
    ######################### Rp2paths and RetroPath2 #####################
    #######################################################################


    @staticmethod
    def rp2ToCollection(rp2_pathways,
                        rp2paths_compounds,
                        rp2paths_pathways,
                        rpcollection_output,
                        upper_flux_bound=999999.0,
                        lower_flux_bound=0.0,
                        max_subpaths_filter=100,
                        pathway_id='rp_pathway',
                        compartment_id='MNXC3',
                        species_group_id='central_species',
                        sink_species_group_id='rp_sink_species',
                        pubchem_search=False,
                        rpcache=None):
        if not rpcache:
            rpcache = rpCache()
            #WARNING: there is no real need to populate the cache as its done dynamically now
            rpcache.populateCache()
        rpreader = rpReader(rpcache)
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            os.mkdir(os.path.join(tmp_output_folder, 'models'))
            #make the log
            rpreader_log = {}
            rpreader_log['rpreader'] = {}
            rpreader_log['rpreader'][time.time()] = {'rp2_pathways': rp2_pathways,
                                                     'rp2paths_compounds': rp2paths_compounds,
                                                     'rp2paths_pathways': rp2paths_pathways,
                                                     'rpcollection_output': rpcollection_output,
                                                     'upper_flux_bound': upper_flux_bound,
                                                     'lower_flux_bound': lower_flux_bound,
                                                     'max_subpaths_filter': max_subpaths_filter,
                                                     'pathway_id': pathway_id,
                                                     'compartment_id': compartment_id,
                                                     'species_group_id': species_group_id,
                                                     'sink_species_group_id': sink_species_group_id,
                                                     'pubchem_search': pubchem_search}
            json.dump(rpreader_log, open(os.path.join(tmp_output_folder, 'log.json'), 'w'))
            #make the run
            status = rpreader.rp2ToSBML(rp2_pathways,
                                        rp2paths_compounds,
                                        rp2paths_pathways,
                                        os.path.join(tmp_output_folder, 'models'),
                                        upper_flux_bound,
                                        lower_flux_bound,
                                        max_subpaths_filter,
                                        pathway_id,
                                        compartment_id,
                                        species_group_id,
                                        sink_species_group_id,
                                        pubchem_search)
            if len(glob.glob(os.path.join(tmp_output_folder, 'models', '*')))==0:
                logging.error('Output has not produced any models')
                return False
            if not rpcollection_output.endswith('.rpcol'):
                rpcollection_output += '.rpcol'
            if os.path.exists(rpcollection_output):
                logging.warning('The path '+str(rpcollection_output)+' already exists... overwriting it')
            #save the whole thing
            with tarfile.open(rpcollection_output, "w:xz") as tar:
                tar.add(tmp_output_folder, arcname='rpsbml_collection')
        return status
    

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
        rp_strc = self.readRp2PathsCompounds(rp2paths_compounds)
        rp_transformations, rp_compounds = self.readRpPathways(rp2_pathways)
        # TODO: make sure that you account for the fact that each reaction may have multiple associated reactions
        rp_paths = self.readRp2PathsPathways(rp2paths_pathways)
        sink_species = []
        # for each line or rp2paths_pathways:
        #     generate comb
        #     for each combinant:
        #         rank
        #         process
        #         add cofactors
        #         dedup
        #### pathToSBML ####
        compartment_id = self._checkCompDeprecated(compartment_id)
        for path_num in rp_paths:
            # first level is the list of lists of sub_steps
            # second is itertools all possible combinations using product
            alt_path_num = 1
            # sopX subpaths of the current rp2path pathway
            list_rpsbml = []
            for comb_path in list(itertools.product(*[[(i,y) for y in rp_paths[path_num][i]] for i in rp_paths[path_num]])):
                steps = []
                for i, y in comb_path:
                    steps.append(rp_paths[path_num][i][y])
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
                path_group = rpsbml.createGroup(pathway_id)
                rpsbml.addUpdateBRSynth(path_group, 'path_id', path_id, None, False, False, False)
                rpsbml.addUpdateBRSynth(path_group, 'sub_path_id', alt_path_num, None, False, False, False)
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
                    #rpsbml = self._addSpecies(rpsbml, cid, rp_compounds, compartment_id, chem_name, spe, species_group_id, sink_species_group_id)
                    #create the species according to the type of group they belong to
                    if not species in self.species_obj_dict:
                        self.species_obj_dict.append(species)
                    if cid in rp_compounds:
                        rpsbml.createSpecies(cid,
                                             compartment_id,
                                             species.name,
                                             species.xref,
                                             species.inchi,
                                             species.inchikey,
                                             species.smiles,
                                             species_group_id,
                                             sink_species_group_id,
                                             use_species_id_as_is=False)
                    else:
                        rpsbml.createSpecies(cid,
                                             compartment_id,
                                             species.name,
                                             species.xref,
                                             species.inchi,
                                             species.inchikey,
                                             species.smiles,
                                             species_group_id,
                                             use_species_id_as_is=False)
                #4) Add the complete reactions and their annotations
                for step in steps:
                    # add the substep to the model
                    #step['sub_step'] = alt_path_num
                    rpsbml.createReaction('RP'+str(step['step']), # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['step']),
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          step,
                                          compartment_id,
                                          rp_transformations[step['transformation_id']]['rule'],
                                          {'ec': rp_transformations[step['transformation_id']]['ec']},
                                          pathway_id,
                                          use_species_id_as_is=False)
                #5) Adding the consumption of the target
                target_step = {'rule_id': None,
                               'left': {[i for i in all_cid if i[:6]=='TARGET'][0]: 1}, #warning this is dangerous
                               'right': {},
                               'step': None,
                               'transformation_id': None,
                               'rule_score': None,
                               'rule_ori_reac': None}
                rpsbml.createReaction('RP1_sink',
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      target_step,
                                      compartment_id,
                                      use_species_id_as_is=False)
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
                self.logger.debug('out_dir: '+str(out_dir))
                rpsbml.writeSBML(out_dir)
        return True


    def readRp2PathsCompounds(self, path):
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
                        res_conv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchi'})
                        rp_strc[row[0]]['inchi'] = res_conv['inchi']
                    except NotImplementedError as e:
                        self.logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                try:
                    if not strc_info:
                        strc_info = self.queryCIDstrc(row[0])
                    rp_strc[row[0]]['inchikey'] = strc_info['inchikey']
                    #try to generate them yourself by converting them directly
                    #TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
                except KeyError:
                    try:
                        res_conv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})
                        rp_strc[row[0]]['inchikey'] = res_conv['inchikey']
                    except NotImplementedError as e:
                        self.logger.warning('Could not convert the following SMILES to InChI key: '+str(row[1]))
        except (TypeError, FileNotFoundError) as e:
            self.logger.error('Could not read the compounds file ('+str(path)+')')
            raise RuntimeError
        return rp_strc


    def readRpPathways(self, path):
        """Function to parse the scope.csv file

        Extract the reaction rules from the retroPath2.0 output using the scope.csv file

        :param path: Path to the compounds file

        :type path: str

        :rtype: tuple
        :return: The RetroPath transformation and the list of sink molecules
        """
        rp_transformations = {}
        rp_compounds = []
        #### we might pass binary in the REST version
        reader = None
        if isinstance(path, bytes):
            reader = csv.reader(StringIO(path.decode('utf-8')), delimiter=',')
        else:
            try:
                reader = csv.reader(open(path, 'r'), delimiter=',')
            except FileNotFoundError:
                self.logger.error('Could not read the compounds file: '+str(path))
                return {}
        next(reader)
        for row in reader:
            if not row[1] in rp_transformations:
                rp_transformations[row[1]] = {}
                rp_transformations[row[1]]['rule'] = row[2]
                rp_transformations[row[1]]['ec'] = [i.replace(' ', '') for i in row[11][1:-1].split(',') if not i.replace(' ', '')=='NOEC']
            if row[7]=='1':
                for i in row[8].replace(']', '').replace('[', '').replace(' ', '').split(','):
                    rp_compounds.append(i)
        #self.logger.info(rp_transformations)
        #self.logger.info(rp_compounds)
        return rp_transformations, list(set(rp_compounds))


    def readRp2PathsPathways(self, rp2paths_pathways):
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
                self.logger.error('Cannot convert path_id to int ('+str(row[0])+')')
                #return {}
                return False
            #################################
            rule_ids = row[2].split(',')
            if rule_ids==None:
                self.logger.warning('The rulesIds is None')
                #pass # or continue
                continue
            ###WARNING: This is the part where we select some rules over others
            # we do it by sorting the list according to their score and taking the topx
            sub_path_step = 1
            for r_id in rule_ids:
                rr_reacts = self.queryRRreactions(r_id)
                if not rr_reacts:
                    continue
                for rea_id in rr_reacts:
                    tmp_reac = {'rule_id': r_id,
                                'rule_ori_reac': rea_id,
                                'rule_score': rr_reacts[rea_id]['rule_score'],
                                'right': {},
                                'left': {},
                                'path_id': int(row[0]),
                                'step': path_step,
                                'transformation_id': row[1][:-2]}
                    ############ LEFT ##############
                    for l in row[3].split(':'):
                        tmp_l = l.split('.')
                        #tmp_reac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                        cid = '' #TODO: change this
                        cid = self._checkCIDdeprecated(tmp_l[1])
                        try:
                            tmp_reac['left'][cid] = int(tmp_l[0])
                        except ValueError:
                            self.logger.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                            #return {}
                            return False
                    ############## RIGHT ###########
                    for r in row[4].split(':'):
                        tmp_r = r.split('.')
                        #tmp_reac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                        cid = '' #TODO change this
                        cid = self._checkCIDdeprecated(tmp_r[1])
                        try:
                            tmp_reac['right'][cid] = int(tmp_r[0])
                        except ValueError:
                            self.logger.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                            return False
                    #################################
                    if not int(row[0]) in rp_paths:
                        rp_paths[int(row[0])] = {}
                    if not int(path_step) in rp_paths[int(row[0])]:
                        rp_paths[int(row[0])][int(path_step)] = {}
                    rp_paths[int(row[0])][int(path_step)][int(sub_path_step)] = tmp_reac
                    sub_path_step += 1
            '''
            tmp_rr_reactions[str(r_id)+'__'+str(rea_id)] = rr_reacts[rea_id]
            for single_rule in tmp_rr_reactions:
                tmp_reac = {'rule_id': single_rule.split('__')[0],
                            'rule_ori_reac': single_rule.split('__')[1],
                            'rule_score': self.rr_reactions[single_rule.split('__')[0]][single_rule.split('__')[1]]['rule_score'],
                            'right': {},
                            'left': {},
                            'path_id': int(row[0]),
                            'step': path_step,
                            'transformation_id': row[1][:-2]}
                ############ LEFT ##############
                for l in row[3].split(':'):
                    tmp_l = l.split('.')
                    #tmp_reac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                    cid = '' #TODO: change this
                    cid = self._checkCIDdeprecated(tmp_l[1])
                    try:
                        tmp_reac['left'][cid] = int(tmp_l[0])
                    except ValueError:
                        self.logger.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                        #return {}
                        return False
                ############## RIGHT ###########
                for r in row[4].split(':'):
                    tmp_r = r.split('.')
                    #tmp_reac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                    cid = '' #TODO change this
                    cid = self._checkCIDdeprecated(tmp_r[1])
                    try:
                        tmp_reac['right'][cid] = int(tmp_r[0])
                    except ValueError:
                        self.logger.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                        return False
                #################################
                if not int(row[0]) in rp_paths:
                    rp_paths[int(row[0])] = {}
                if not int(path_step) in rp_paths[int(row[0])]:
                    rp_paths[int(row[0])][int(path_step)] = {}
                rp_paths[int(row[0])][int(path_step)][int(sub_path_step)] = tmp_reac
                sub_path_step += 1
            '''
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
        if not species.cid and species.inchikey:
            species.cid = self.queryInchiKeyCID(species.inchikey)
        if not species.cid:
            self.logger.error('There must be a cid of the Species object')
            return False
        #common name of the species
        if not species.name:
            species.name = self.queryCIDname(species.cid)
        ###Try to fill from cache
        if not species.xref:
            species.xref = self.queryCIDxref(species.cid)
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
                    species.xref = self.queryCIDxref(self.queryChebiCID(pubres['xref']['chebi'][0]))
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


    #############################################################################################
    ############################### COFACTORS ###################################################
    #############################################################################################


    def _completeReac(self, step, rr_reac, full_reac, mono_side, rr_string, pathway_cmp):
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
                self.logger.warning('There should be only one compound on the left for monocomponent reaction: rr_mono_cmp: '+str(rr_mono_cmp)+' step_mono_cmp: '+str(step_mono_cmp))
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
        for to_add in full_reac.keys()-rr_reac.keys():
            step.update({to_add: full_reac[to_add]})
            ### update the reaction rule string
            try:
                smi = self.queryCIDstrc(to_add)['smiles']
                if not smi==None:
                    for sto_add in range(int(full_reac[to_add])):
                        rr_string += '.'+str(smi)
            except KeyError:
                self.logger.warning('Cannot find smiles structure for '+str(to_add))
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
                        self.logger.warning('full_reac stochio should never be smaller than step')
                        continue
                    for i in range(stochio_diff):
                        ### update the reaction rule string
                        try:
                            smi = self.queryCIDstrc(step_spe)['smiles']
                            if not smi==None:
                                rr_string += '.'+str(smi)
                        except KeyError:
                            self.logger.warning('Cannot find smiles structure for '+str(step_spe))
            elif step_spe in pathway_cmp:
                if pathway_cmp[step_spe] in full_reac:
                    if not step[step_spe]==full_reac[pathway_cmp[step_spe]]:
                        step[step_spe] = full_reac[pathway_cmp[step_spe]]
            #Its fine if the stochio is not updated, better than ignoring a whole pathway
                #else:
                #    self.logger.warning('Cannot find '+str(step_spe)+' in full reaction')
                #    return False
            #else:
            #    self.logger.warning('Cannot find '+str(step_spe)+' in pathway_cmp')
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
                isSuccess, reac_smiles_left = self._completeReac(step['right'],
                                                                 self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                                                                 self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'])]['right'],
                                                                 True,
                                                                 reac_smiles_left,
                                                                 pathway_cmp)
                if not isSuccess:
                    self.logger.warning('Could not recognise reaction rule for step (1): '+str(step))
                    return False
            except KeyError:
                self.logger.warning('Could not find the full reaction for reaction (1): '+str(step))
                return False
            try:
                isSuccess, reac_smiles_right = self._completeReac(step['left'],
                                                                  self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                                                                  self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'])]['left'],
                                                                  False,
                                                                  reac_smiles_right,
                                                                  pathway_cmp)
                if not isSuccess:
                    self.logger.warning('Could not recognise reaction rule for step (2): '+str(step))
                    return False
            except KeyError:
                self.logger.warning('Could not find the full reaction for reaction (2): '+str(step))
                return False
        elif self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==1:
            try:
                isSuccess, reac_smiles_left = self._completeReac(step['right'],
                                                                 self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                                                                 self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'])]['left'],
                                                                 True,
                                                                 reac_smiles_left,
                                                                 pathway_cmp)
                if not isSuccess:
                    self.logger.error('Could not recognise reaction rule for step (3): '+str(step))
                    return False
            except KeyError:
                self.logger.warning('Could not find the full reaction for reaction (3): '+str(step))
                return False
            try:
                isSuccess, reac_smiles_right = self._completeReac(step['left'],
                                                                  self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                                                                  self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'])]['right'],
                                                                  False,
                                                                  reac_smiles_right,
                                                                  pathway_cmp)
                if not isSuccess:
                    self.logger.error('Could not recognise reaction rule for step (4): '+str(step))
                    return False
            except KeyError:
                self.logger.warning('Could not find the full reaction for reaction (4): '+str(step))
                return False
        else:
            self.logger.error('Relative direction can only be 1 or -1: '+str(self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']))
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
        rpsbml_dict = rpsbml.asDict(pathway_id)
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
                                             species.smiles,
                                             use_species_id_as_is=False)
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
                self.logger.warning('Cannot find cofactors... skipping')
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
    # @param in_file The input JSON file
    # @param mnx_header Reorganise the results around the target MNX products
    # @return Dictionnary of SBML
    def _parseTSV(self, in_file, remove_inchi_4p=False, mnx_header=False):
        data = {}
        try:
            for row in csv.DictReader(open(in_file), delimiter='\t'):
                ######## path_id ######
                try:
                    path_id = int(row['pathway_ID'])
                except ValueError:
                    self.logger.error('Cannot convert pathway ID: '+str(row['pathway_ID']))
                    continue
                if not path_id in data:
                    data[path_id] = {}
                    data[path_id]['is_valid'] = True
                    data[path_id]['steps'] = {}
                ####### target #########
                if not 'target' in data[path_id]:
                    data[path_id]['target'] = {}
                    data[path_id]['target']['name'] = row['target_name']
                    if remove_inchi_4p:
                        data[path_id]['target']['inchi'] = '/'.join([row['target_structure'].split('/')[i] for i in range(len(row['target_structure'].split('/'))) if i<4])
                    else:
                        data[path_id]['target']['inchi'] = row['target_structure']
                ####### step #########
                try:
                    step_id = int(row['step'])
                except ValueError:
                    self.logger.error('Cannot convert step ID: '+str(row['step']))
                    data[path_id]['is_valid'] = False
                    continue
                if step_id==0:
                    continue
                elif step_id==1:
                    data[path_id]['organism'] = row['organism'].replace(' ', '')
                    data[path_id]['reference'] = row['reference'].replace(' ', '')
                data[path_id]['steps'][step_id] = {}
                ##### substrates #########
                data[path_id]['steps'][step_id]['substrates'] = []
                len_db_ref = len(row['substrate_dbref'].split(';'))
                for i in row['substrate_dbref'].split(';'):
                    if i=='':
                        len_db_ref -= 1
                len_strc = len(row['substrate_structure'].split('_'))
                for i in row['substrate_structure'].split('_'):
                    if i=='':
                        len_strc -= 1
                len_sub = len(row['substrate_name'].split(';'))
                for i in row['substrate_name'].split(';'):
                    if i=='':
                        len_sub -= 1
                if len_sub==len_strc==len_sub:
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
                                self.logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                                data[path_id]['is_valid'] = False
                        data[path_id]['steps'][step_id]['substrates'].append(tmp)
                else:
                    self.logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['substrate_name'])+' <--> '+str(row['substrate_structure'])+' <--> '+str(row['substrate_dbref']))
                    data[path_id]['is_valid'] = False
                    continue
                ##### products #########
                data[path_id]['steps'][step_id]['products'] = []
                len_db_ref = len(row['product_dbref'].split(';'))
                for i in row['product_dbref'].split(';'):
                    if i=='':
                        len_db_ref -= 1
                len_strc = len(row['product_structure'].split('_'))
                for i in row['product_structure'].split('_'):
                    if i=='':
                        len_strc -= 1
                len_sub = len(row['product_name'].split(';'))
                for i in row['product_name'].split(';'):
                    if i=='':
                        len_sub -= 1
                if len_sub==len_strc==len_db_ref:
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
                                data[path_id]['is_valid'] = False
                                self.logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                        data[path_id]['steps'][step_id]['products'].append(tmp)
                else:
                    self.logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['product_name'])+' <--> '+str(row['product_structure'])+' <--> '+str(row['product_dbref']))
                    data[path_id]['is_valid'] = False
                if not row['uniprot']=='':
                    data[path_id]['steps'][step_id]['uniprot'] = row['uniprot'].replace(' ', '').split(';')
                if not row['EC_number']=='':
                    data[path_id]['steps'][step_id]['ec_numbers'] = [i.replace(' ', '') for i in row['EC_number'].split(';')]
                data[path_id]['steps'][step_id]['enzyme_id'] = [i.replace(' ', '') for i in row['enzyme_identifier'].split(';')]
                data[path_id]['steps'][step_id]['enzyme_name'] = row['enzyme_name'].split(';')
        except FileNotFoundError:
            self.logger.error('Cannot open the file: '+str(in_file))
        #now loop through all of them and remove the invalid paths
        to_ret = copy.deepcopy(data)
        for path_id in data.keys():
            if to_ret[path_id]['is_valid']==False:
                del to_ret[path_id]
            else:
                del to_ret[path_id]['is_valid']
        #reorganise the results around the target products mnx
        if not mnx_header:
            return to_ret
        else:
            to_ret_two = {}
            for path_id in to_ret:
                try:
                    final_pro_mnx = to_ret[path_id]['steps'][max(to_ret[path_id]['steps'])]['products'][0]['dbref']['mnx'][0]
                except KeyError:
                    self.logger.error('The species '+str(to_ret[path_id]['steps'][max(to_ret[path_id]['steps'])]['products'][0]['name'])+' does not contain a mnx database reference... skipping whole pathway number '+str(path_id))
                    #continue
                if not final_pro_mnx in to_ret_two:
                    to_ret_two[final_pro_mnx] = {}
                to_ret_two[final_pro_mnx][path_id] = to_ret[path_id]
            return to_ret_two


    ## Parse the validation TSV to SBML
    #
    # Parse the TSV file to SBML format and adds them to the sbml_paths
    #
    # @param self Object pointer
    # @param in_file Input file
    # @param compartment_id compartment of the
    # TODO: update this with the new SBML groups
    def TSVtoSBML(self,
                  in_file,
                  output_folder=None,
                  upper_flux_bound=99999.0,
                  lower_flux_bound=0.0,
                  compartment_id='MNXC3',
                  pathway_id='rp_pathway',
                  species_group_id='central_species',
                  header_name=''):
        data = self._parseTSV(in_file)
        sbml_paths = {}
        if header_name=='':
            header_name = in_file.split('/')[-1].replace('.tsv', '').replace('.csv', '')
        # TODO: need to exit at this loop
        for path_id in data:
            try:
                compartment_id = self._checkCompDeprecated(compartment_id)
            except KeyError:
                self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
                return False
            rpsbml = rpSBML.rpSBML(model_name=header_name+'_'+str(path_id), rpcache=self.rpcache)
            # 1) create a generic Model, ie the structure and unit definitions that we will use the most
            ##### TODO: give the user more control over a generic model creation:
            # -> special attention to the compartment
            rpsbml.genericModel(header_name+'_Path'+str(path_id),
                                header_name+'_Path'+str(path_id),
                                self.queryCompXref(compartment_id),
                                compartment_id,
                                upper_flux_bound,
                                lower_flux_bound)
            # 2) create the pathway (groups)
            rpsbml.createGroup(pathway_id)
            rpsbml.createGroup(species_group_id)
            # 3) find all the unique species and add them to the model
            all_chem = []
            for step_num in data[path_id]['steps']:
                # because of the nature of the input we need to remove duplicates
                for i in data[path_id]['steps'][step_num]['substrates']+data[path_id]['steps'][step_num]['products']:
                    if not i in all_chem:
                        all_chem.append(i)
            # add them to the SBML
            for chem in all_chem:
                # PROBLEM: as it stands one expects the cid to be MNX
                if 'mnx' in chem['dbref']:
                    # must list the different models
                    cid = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                else:
                    # TODO: add the species with other types of xref in annotation
                    self.logger.warning('Some species are not referenced by a MNX id and will be ignored')
                    # try CHEBI
                    try:
                        cid = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                        cid = 'CHEBI_'+str(cid)
                    except KeyError:
                        # TODO: need to find a better way
                        self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                        tmp_db_name = list(chem['dbref'].keys())[0]
                        cid = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                        cid = str(tmp_db_name)+'_'+str(cid)
                    # break
                # try to conver the inchi into the other structures
                smiles = None
                inchikey = None
                try:
                    res_conv = self._convert_depiction(idepic=chem['inchi'], itype='inchi', otype={'smiles','inchikey'})
                    smiles = res_conv['smiles']
                    inchikey = res_conv['inchikey']
                except NotImplementedError as e:
                    self.logger.warning('Could not convert the following InChI: '+str(chem['inchi']))
                # create a new species
                try:
                    chem_name = self.queryCIDname(cid)['name']
                except KeyError:
                    chem_name = cid
                # compile as much info as you can
                # xref
                try:
                    # TODO: add the xref from the document
                    spe_xref = self.queryCIDxref(cid)
                except KeyError:
                    #spe_xref = {}
                    spe_xref = chem['dbref']
                # inchi
                strc_info = None
                try:
                    if not strc_info:
                        strc_info = self.queryCIDstrc(cid)
                    spe_inchi = strc_info['inchi']
                except KeyError:
                    spe_inchi = chem['inchi']
                # inchikey
                try:
                    if not strc_info:
                        strc_info = self.queryCIDstrc(cid)
                    spe_inchikey = strc_info['inchikey']
                except KeyError:
                    spe_inchikey =  res_conv['inchikey']
                # smiles
                try:
                    if not strc_info:
                        strc_info = self.queryCIDstrc(cid)
                    spe_smiles = strc_info['smiles']
                except KeyError:
                    spe_smiles = res_conv['smiles']
                # pass the information to create the species
                rpsbml.createSpecies(cid,
                                     compartment_id,
                                     chem_name,
                                     spe_xref,
                                     spe_inchi,
                                     spe_inchikey,
                                     spe_smiles,
                                     species_group_id,
                                     use_species_id_as_is=False)
            # 4) add the complete reactions and their annotations
            # create a new group for the measured pathway
            # need to convert the validation to step for reactions
            for step_num in data[path_id]['steps']:
                to_send = {'left': {}, 'right': {}, 'rule_id': None, 'rule_ori_reac': None, 'rule_score': None, 'path_id': path_id, 'step': step_num, 'sub_step': None}
                for chem in data[path_id]['steps'][step_num]['substrates']:
                    if 'mnx' in chem['dbref']:
                        cid = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        # try CHEBI
                    else:
                        self.logger.warning('Not all the species to have a MNX ID')
                        # break
                        try:
                            cid = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            cid = 'CHEBI_'+str(cid)
                        except KeyError:
                            # TODO: need to find a better way
                            self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmp_db_name = list(chem['dbref'].keys())[0]
                            cid = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            cid = str(tmp_db_name)+'_'+str(cid)
                    to_send['left'][cid] = 1
                for chem in data[path_id]['steps'][step_num]['products']:
                    if 'mnx' in chem['dbref']:
                        cid = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        # try CHEBI
                    else:
                        self.logger.warning('Need all the species to have a MNX ID')
                        try:
                            cid = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            cid = 'CHEBI_'+str(cid)
                        except KeyError:
                            # TODO: need to find a better way
                            self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmp_db_name = list(chem['dbref'].keys())[0]
                            cid = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            cid = str(tmp_db_name)+'_'+str(cid)
                    to_send['right'][cid] = 1
                        # break
                # if all are full add it
                reac_xref = {}
                if 'ec_numbers' in data[path_id]['steps'][step_num]:
                    reac_xref['ec'] = data[path_id]['steps'][step_num]['ec_numbers']
                if 'uniprot' in data[path_id]['steps'][step_num]:
                    reac_xref['uniprot'] = data[path_id]['steps'][step_num]['uniprot']
                self.logger.debug('#########################################')
                self.logger.debug(to_send)
                self.logger.debug('#########################################')
                rpsbml.createReaction(header_name+'_Step'+str(step_num),
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      to_send,
                                      compartment_id,
                                      None,
                                      reac_xref,
                                      pathway_id,
                                      use_species_id_as_is=False)
                if step_num==1:
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
                    for chem in data[path_id]['steps'][step_num]['products']:
                        try:
                            # smallest MNX
                            cid = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        except KeyError:
                            # try CHEBI
                            try:
                                cid = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                                cid = 'CHEBI_'+str(cid)
                            except KeyError:
                                self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                                tmp_db_name = list(chem['dbref'].keys())[0]
                                cid = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                                cid = str(tmp_db_name)+'_'+str(cid)
                        target_step['left'][cid] = 1
                    rpsbml.createReaction(header_name+'_Step1_sink',
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          target_step,
                                          compartment_id,
                                          use_species_id_as_is=False)
                    rpsbml.createFluxObj('rpFBA_obj', header_name+'_Step1_sink', 1, True)
            if output_folder:
                rpsbml.writeSBML(output_folder)
            else:
                sbml_paths[header_name+'_Path'+str(path_id)] = rpsbml
        if output_folder:
            return {}
        else:
            return sbml_paths
