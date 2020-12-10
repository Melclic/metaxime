##############################################################
########################## PIPELINE ##########################
##############################################################

from equilibrator_api import ComponentContribution

from metaxime import rpCache
from metaxime import rpReader
from metaxime import rpFBA
from metaxime import rpEquilibrator
from metaxime import rpSelenzyme
from metaxime import rpGlobalScore
from metaxime import rpExtractSink

import callRP2
import callRP2paths
import callRR

from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs

model_list = {'b_subtilis_iYO844': '/home/models/b_subtilis_iYO844.sbml',
              'e_coli_iJO1366': '/home/models/e_coli_iJO1366.sbml',
              'p_putida_iJN746': '/home/models/p_putida_iJN746.sbml',
              'e_coli_core_model': '/home/models/e_coli_core_model.sbml',
              'e_coli_iJR904': '/home/models/e_coli_iJR904.sbml',
              's_cerevisiae_iMM904': '/home/models/s_cerevisiae_iMM904.sbml',
              'y_lipolytica_iMK735': '/home/models/y_lipolytica_iMK735.sbml',
              'e_coli_iAF1260': '/home/models/e_coli_iAF1260.sbml',
              'e_coli_iML1515': '/home/models/e_coli_iML1515.sbml',
              's_cerevisiae_iND750': '/home/models/s_cerevisiae_iND750.sbml'}

model_taxo_id = {'b_subtilis_iYO844': 1423,
                 'e_coli_iJO1366': 83333,
                 'p_putida_iJN746': 160488,
                 'e_coli_core_model': 83333,
                 'e_coli_iJR904': 83333,
                 's_cerevisiae_iMM904': 559292,
                 'y_lipolytica_iMK735': 284591,
                 'e_coli_iAF1260': 83333,
                 'e_coli_iML1515': 83333,
                 's_cerevisiae_iND750': 559292}


global_rpcache = rpCache()
global_rpcache.populateCache()
global_cc = ComponentContribution()
#global_selenzyme = rpSelenzyme(cache_tar_path=os.path.dirname('home', 'metaxime', 'metaxime', 'input_cache', 'rpselenzyme_data.tar.xz'))


def pipeline(rpcollection_file,
             target_smiles,
             gem_name,
             max_steps,
             rules_diameters='2,4,6,8,10,12,14,16',
             rules_type='all',
             topx=100,
             timeout=90.0,
             partial_retro=False,
             ph=7.5,
             ionic_strength=200,
             temp_k=298.15):
    target_inchi = MolToInchi(MolFromSmiles(target_smiles, sanitize=True))
    with tempfile.TemporaryDirectory() as tmp_output_folder:
        ############# source file #############
        source_file = os.path.join(tmp_output_folder, 'source.csv')
        with open(source_file, 'w') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(['Name', 'InChI'])
            filewriter.writerow(['target', target_inchi])
        ############ sink file ##############
        try:
            gem_file = model_list[gem_name]
            taxo_id = model_taxo_id[gem_name]
        except KeyError:
            logging.error('Cannot find the following GEM model: '+str(gem_name))
            return False, 'gem'
        sink_file = os.path.join(tmp_output_folder, 'sink.csv')
        gensink_status = rpExtractSink.genSink(gem_file, sink_file, remove_dead_end=True)
        if not gensink_status:
            return False, 'gensink'
        ########## reaction rules ##########
        rules_file = os.path.join(tmp_output_folder, 'reaction_rules.csv')
        rr_status = callRR.passRules(rules_file, rules_type, rules_diameters, 'csv')
        if not rr_status:
            return False, 'reactionrules'
        ########## Retropath2 ##############
        rp2_file = os.path.join(tmp_output_folder, 'rp2_pathways.csv')
        rp2_status = callRP2.run(rp2_file,
                                 source_file,
                                 sink_file,
                                 rules_file,
                                 max_steps,
                                 topx=topx,
                                 timeout=timeout,
                                 partial_retro=partial_retro)
        if not rp2_status:
            return False, 'rp2'
        ######### rp2paths ################
        rp2paths_pathways_file = os.path.join(tmp_output_folder, 'rp2paths_pathways.csv')
        rp2paths_compounds_file = os.path.join(tmp_output_folder, 'rp2paths_compounds.tsv')
        rp2paths_status = callRP2paths.run(rp2_file, rp2paths_pathways_file, rp2paths_compounds_file)
        if not rp2paths_status:
            return False, 'rp2paths'
        ####### Analysis #################
        rpreader_status = rpReader.rp2ToCollection(rp2_file,
                                                   rp2paths_compounds_file,
                                                   rp2paths_pathways_file,
                                                   rpcollection_file,
                                                   rpcache=global_rpcache)
        if not rpreader_status:
            return False, 'rpreader'
        rpfba_status = rpFBA.runCollection(rpcollection_file,
                                           gem_file,
                                           rpcollection_file,
                                           num_workers=1,
                                           keep_merged=True,
                                           del_sp_pro=False,
                                           del_sp_react=False,
                                           rpcache=global_rpcache)
        if not rpfba_status:
            return False, 'rpfba'
        rpeq_status = rpEquilibrator.runCollection(rpcollection_file,
                                                   rpcollection_file,
                                                   cc=global_cc,
                                                   ph=ph,
                                                   ionic_strength=ionic_strength,
                                                   temp_k=temp_k,
                                                   rpcache=global_rpcache)
        if not rpeq_status:
            return False, 'rpeq'
        rpsel_status = rpSelenzyme.runCollection(rpcollection_file,
                                                 taxo_id,
                                                 rpcollection_file,
                                                 cache_path=os.path.join('/home', 'metaxime', 'input_cache', 'rpselenzyme_data.tar.xz'),
                                                 rpcache=global_rpcache)
        if not rpsel_status:
            return False, 'rpsel'
        rpglo_status = rpGlobalScore.runCollection(rpcollection_file,
                                                   rpcollection_file,
                                                   rpcache=global_rpcache)
        if not rpglo_status:
            return False, 'rpglo'
        return True, 'success'

