import json
import logging
import os
import glob

from metaxime import rpCache
from metaxime import rpGraph
from metaxime import rpReader
from metaxime import rpFBA
from metaxime import rpSBML
from metaxime import rpEquilibrator
from metaxime import rpSelenzyme
from metaxime import rpGlobalScore

#TODO: delete this file.... temporary

logging.basicConfig(
    level=logging.DEBUG,
    #level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


model_taxo_id = {'b_subtilis_iYO844': 224308,
                 'e_coli_iJO1366': 83333,
                 'p_putida_iJN746': 160488,
                 'e_coli_core_model': 83333,
                 'e_coli_iJR904': 83333,
                 's_cerevisiae_iMM904': 559292,
                 'e_coli_iAF1260': 83333,
                 'y_lipolytica_iMK735': 284591,
                 's_cerevisiae_iND750': 559292,
                 'e_coli_iML1515': 83333}

global_rpcache = rpCache()
global_rpcache.populateCache()

ph = 7.5
ionic_strength = 200
temp_k = 298.15

selenzyme_cache_file = '/Users/melchior/workspace/melclic/metaxime/metaxime/input_cache/rpselenzyme_data.tar.xz'


for res_dir in glob.glob('/Users/melchior/workspace/melclic/metaxime/laser_compounds_run/*'):
    for strain_dir in glob.glob(os.path.join(res_dir, '*')):
        logging.info('####### '+str(strain_dir)+' #######')
        rp2_file = os.path.join(strain_dir, 'rp_pathways.csv')
        rp2paths_compounds_file = os.path.join(strain_dir, 'rp2paths_compounds.csv')
        rp2paths_pathways_file = os.path.join(strain_dir, 'rp2paths_pathways.csv')
        rpcollection_file = os.path.join(strain_dir, 'rp_collection.tar.xz')
        gem_file = os.path.join(strain_dir, 'model.sbml')
        if not os.path.exists(rp2_file) or not os.path.exists(rp2paths_compounds_file) or not os.path.exists(rp2paths_pathways_file):
            logging.warning('Either RP2 or RP2paths was not succesfull')
            logging.warning('rp2_file: '+str(rp2_file))
            logging.warning('rp2paths_pathways_file: '+str(rp2paths_pathways_file))
            logging.warning('rp2paths_compounds_file: '+str(rp2paths_compounds_file))
            continue
        ##################################
        ####### Analysis #################
        ##################################
        logging.info('---------- rpReader -----')
        rpreader_status = rpReader.rp2ToCollection(rp2_file,
                                                   rp2paths_compounds_file,
                                                   rp2paths_pathways_file,
                                                   rpcollection_file,
                                                   rpcache=global_rpcache)
        if not rpreader_status:
            logging.error('Problem running rpReader')
            continue
        logging.info('---------- rpFBA -----')
        if not os.path.exists(gem_file):
            logging.error('Cannot find the GEM model: '+str(gem_file))
            continue
        rpfba_status = rpFBA.runCollection(rpcollection_file,
                                           gem_file,
                                           rpcollection_file,
                                           num_workers=1,
                                           keep_merged=True,
                                           del_sp_pro=False,
                                           del_sp_react=False,
                                           rpcache=global_rpcache)
        if not rpfba_status:
            logging.error('Problem running rpFBA')
            continue
        logging.info('---------- rpEquilibrator -----')
        rpeq_status = rpEquilibrator.runCollection(rpcollection_file,
                                                   rpcollection_file,
                                                   ph=ph,
                                                   ionic_strength=ionic_strength,
                                                   temp_k=temp_k,
                                                   rpcache=global_rpcache)
        if not rpeq_status:
            logging.error('Problem running rpEquilibrator')
            continue
        logging.info('---------- rpSelenzyme -----')
        #### try to revocer the tax_id
        taxo_id = None
        try:
            with open(os.path.join(strain_dir, 'rp2_report.txt'), 'r') as rp2_report:
                lines = []
                for i in rp2_report:
                    lines.append(i)
                s = lines[2].replace('GEM SBML: ', '').replace('.sbml\n', '')
                taxo_id = model_taxo_id[s]
        except IndexError:
            logging.error('Index error in lines: '+str(lines))
            continue
        except FileNotFoundError:
            logging.error('Cannot find the file')
            continue
        except KeyError:
            logging.error('KeyError for '+str(s))
            continue
        if not taxo_id:
            logging.error('Taxonomy id has not been defined')
            continue
        rpsel_status = rpSelenzyme.runCollection(rpcollection_file,
                                                 taxo_id,
                                                 rpcollection_file,
                                                 cache_path=selenzyme_cache_file,
                                                 rpcache=global_rpcache)
        if not rpsel_status:
            logging.error('Problem running rpSelenzyme')
            continue
        logging.info('---------- rpGlobalScore -----')
        rpglo_status = rpGlobalScore.runCollection(rpcollection_file,
                                                   rpcollection_file,
                                                   rpcache=global_rpcache)
        if not rpglo_status:
            logging.error('Problem running rpGlobalScore')
            continue
        ################ generate the various netowrk JSON and SBML JSON #####
        logging.info('------------ JSON ----------')
        rpSBML.batchAsDict(rpcollection_file)
        rpGraph.batchNetworkDict(rpcollection_file, is_gem_sbml=False)
        logging.info('-------- done --------')
