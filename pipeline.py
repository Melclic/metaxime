##############################################################
########################## PIPELINE ##########################
##############################################################

import os
import csv
import subprocess
import glob
import sys
import resource
import tempfile
import json
import sys
import argparse
import tarfile
import shutil
import os
import rq
from equilibrator_api import ComponentContribution


from metaxime import rpCache
from metaxime import rpGraph
from metaxime import rpReader
from metaxime import rpFBA
from metaxime import rpSBML
from metaxime import rpEquilibrator
from metaxime import rpSelenzyme
from metaxime import rpGlobalScore
from metaxime import rpExtractSink

#import callRP2
#import callRP2paths
#import callRR

from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs

#make sure that you call the logging after metaxime to overwite their configuation
import logging
import logging.config
from logsetup import LOGGING_CONFIG

logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)


KPATH = '/usr/local/knime/knime'
RP_WORK_PATH = '/home/rp2/RetroPath2.0.knwf'
MAX_VIRTUAL_MEMORY = 30000*1024*1024 #30GB
#MAX_VIRTUAL_MEMORY = 30*1024*1024 # 30 GB -- define what is the best
#MAX_VIRTUAL_MEMORY = 3.0*1024*1024 # 1 GB -- define what is the best


model_list = {'b_subtilis_iYO844': '/home/models/b_subtilis_iYO844.sbml',
              'e_coli_iJO1366': '/home/models/e_coli_iJO1366.sbml',
              'p_putida_iJN746': '/home/models/p_putida_iJN746.sbml',
              'e_coli_core_model': '/home/models/e_coli_core_model.sbml',
              'e_coli_iJR904': '/home/models/e_coli_iJR904.sbml',
              's_cerevisiae_iMM904': '/home/models/s_cerevisiae_iMM904.sbml',
              'e_coli_iAF1260': '/home/models/e_coli_iAF1260.sbml',
              'y_lipolytica_iMK735': '/home/models/y_lipolytica_iMK735.sbml',
              's_cerevisiae_iND750': '/home/models/s_cerevisiae_iND750.sbml',
              'e_coli_iML1515': '/home/models/e_coli_iML1515.sbml'}

sink_list = {'b_subtilis_iYO844': '/home/sinks/b_subtilis_iYO844__sink.csv',
              'e_coli_iJO1366': '/home/sinks/e_coli_iJO1366__sink.csv',
              'p_putida_iJN746': '/home/sinks/p_putida_iJN746__sink.csv',
              'e_coli_core_model': '/home/sinks/e_coli_core_model__sink.csv',
              'e_coli_iJR904': '/home/sinks/e_coli_iJR904__sink.csv',
              's_cerevisiae_iMM904': '/home/sinks/s_cerevisiae_iMM904__sink.csv',
              'e_coli_iAF1260': '/home/sinks/e_coli_iAF1260__sink.csv',
              'y_lipolytica_iMK735': '/home/sinks/y_lipolytica_iMK735__sink.csv',
              's_cerevisiae_iND750': '/home/sinks/s_cerevisiae_iND750__sink.csv',
              'e_coli_iML1515': '/home/sinks/e_coli_iML1515__sink.csv'}

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

#global_selenzyme = rpSelenzyme(cache_tar_path=os.path.dirname('home', 'metaxime', 'metaxime', 'input_cache', 'rpselenzyme_data.tar.xz'))

'''
logging.basicConfig(
    level=logging.DEBUG,
    #level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)
'''

def modResJSON(job_id,
               job_status=None,
               job_meta=None,
               input_name=None,
               input_smiles=None,
               input_gem=None,
               input_steps=None,
               output_path=None,
               output_tar=None):
    ####### create/append json #####
    if not os.path.exists('/mx-results/job_results.json'):
        with open('/mx-results/job_results.json', 'w') as jr:
            json.dump({}, jr)
    mx_res = None
    with open('/mx-results/job_results.json', 'r') as jr:
        mx_res = json.load(jr)
    id_pos = {}
    for i in range(len(mx_res)):
        id_pos[mx_res[i]['job']['id']] = i
    if not job_id in id_pos:
        logger.warning('First time creating: '+str(job_id))
        mx_res.append({'job': {'status': job_status, 'meta': job_meta, 'id': job_id},
                       'input': {'name': input_name, 'smiles': input_smiles, 'gem': input_gem, 'steps': input_steps},
                       'output': {'path': output_path, 'tar': output_tar}})
    else:
        if job_status:
            mx_res[id_pos[job_id]]['job']['status'] = job_status
        if job_meta:
            mx_res[id_pos[job_id]]['job']['meta'] = job_meta
        if input_name:
            mx_res[id_pos[job_id]]['input']['name'] = input_name
        if input_smiles:
            mx_res[id_pos[job_id]]['input']['smiles'] = input_smiles
        if input_gem:
            mx_res[id_pos[job_id]]['input']['gem'] = input_gem
        if input_steps:
            mx_res[id_pos[job_id]]['input']['steps'] = input_steps
        if output_path:
            mx_res[id_pos[job_id]]['output']['path'] = output_path
        if output_tar:
            mx_res[id_pos[job_id]]['output']['tar'] = output_tar
    with open('/mx-results/job_results.json', 'w') as jr:
        json.dump(mx_res, jr)



def limit_virtual_memory():
    """Limit the virtual of the subprocess call
    """
    resource.setrlimit(resource.RLIMIT_AS, (MAX_VIRTUAL_MEMORY, resource.RLIM_INFINITY))


def pipeline(target_smiles,
             gem_name,
             max_steps,
             rules_diameters='2,4,6,8,10,12,14,16',
             rules_type='all',
             topx=100,
             sub_timeout=90.0,
             partial_retro=False,
             ph=7.5,
             ionic_strength=200,
             temp_k=298.15):
    logger.info('####### pipieline ######')
    job = rq.get_current_job()
    job.meta['progress'] = 1
    job.save_meta()
    modResJSON(job.id, job_meta=job.meta, job_status='running')
    logger.debug('cache')
    global_rpcache = rpCache()
    logger.debug('populate cache')
    global_rpcache.populateCache()
    logger.info('target_smiles: '+str(target_smiles))
    logger.info('gem_name: '+str(gem_name))
    logger.info('max_steps: '+str(max_steps))
    logger.info('job: '+str(job))
    logger.info('job.id: '+str(job.id))
    target_inchi = MolToInchi(MolFromSmiles(target_smiles, sanitize=True))
    logger.info('target_inchi: '+str(target_inchi))
    with tempfile.TemporaryDirectory() as tmp_output_folder:
        try:
            ############# source file #############
            logger.info('------ source file -----')
            job.meta['progress'] = 2
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            rpcollection_file = os.path.join(tmp_output_folder, 'rpcollection.tar.xz')
            rp_results_folder = os.path.join(tmp_output_folder, 'rp_results')
            os.mkdir(rp_results_folder)
            source_file = os.path.join(rp_results_folder, 'source.csv')
            with open(source_file, 'w') as csvfile:
                filewriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                filewriter.writerow(['Name', 'InChI'])
                filewriter.writerow(['target', target_inchi])
            ############ sink file ##############
            logger.debug('-------- sink file --------')
            job.meta['progress'] = 3
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            try:
                gem_file = model_list[gem_name]
                taxo_id = model_taxo_id[gem_name]
                sink_file = sink_list[gem_name]
            except KeyError:
                logger.error('Cannot find the following GEM model: '+str(gem_name))
                job.meta['err_msg'] = 'gem'
                job.save_meta()
                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                #return False, 'gem', b''
                return False, 'gem', ''
            '''#only when passing an SBML file to process
            sink_file = os.path.join(tmp_output_folder, 'sink.csv')
            gensink_status = rpExtractSink.genSink(gem_file, sink_file, remove_dead_end=True)
            if not gensink_status:
                logger.error('Problem generating sink')
                return False, 'gensink'
            '''
            shutil.copy2(gem_file, os.path.join(rp_results_folder, 'gem_sbml.xml'))
            shutil.copy2(sink_file, os.path.join(rp_results_folder, 'sink.csv'))
            ####################################
            ########## Reaction Rules ##########
            ####################################
            logger.info('---------- reaction rules -----')
            job.meta['progress'] = 4
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            rules_file = os.path.join(rp_results_folder, 'reaction_rules.csv')
            rule_file = None
            if rules_type=='all':
                rule_file = '/home/retrorules/rules_rall_rp2.csv'
            elif rules_type=='forward':
                rule_file = '/home/retrorules/rules_rall_rp2_forward.csv'
            elif rules_type=='retro':
                rule_file = '/home/retrorules/rules_rall_rp2_retro.csv'
            else:
                logger.error('RR: Cannot detect input: '+str(rules_type))
                job.meta['err_msg'] = 'rr_inputerror'
                job.save_meta()
                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                #return False, 'rr_inputerror', b''
                return False, 'rr_inputerror', ''
            #check the input diameters are valid #
            try:
                s_diameters = [int(i) for i in rules_diameters.split(',')]
                valid_diameters = []
                for i in s_diameters:
                    if i not in [2,4,6,8,10,12,14,16]:
                        logger.warning('RR: Diameters must be either 2,4,6,8,10,12,14,16. Ignoring entry: '+str(i))
                    else:
                        valid_diameters.append(i)
            except ValueError:
                logger.error('RR: Invalid diamter entry. Must be int of either 2,4,6,8,10,12,14,16')
                job.meta['err_msg'] = 'rr_invaliddiameters'
                job.save_meta()
                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                #return False, 'rr_invaliddiameters', b''
                return False, 'rr_invaliddiameters', ''
            ##### create temp file to write ####
            with tempfile.TemporaryDirectory() as tmp_rr_folder:
                outfile_path = os.path.join(tmp_rr_folder, 'tmp_rules.csv')
                with open(rule_file, 'r') as rf:
                    with open(outfile_path, 'w') as o:
                        rf_csv = csv.reader(rf)
                        o_csv = csv.writer(o, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
                        o_csv.writerow(next(rf_csv))
                        for row in rf_csv:
                            try:
                                if int(row[4]) in valid_diameters:
                                    o_csv.writerow(row)
                            except ValueError:
                                logger.error('RR: Cannot convert diameter to integer: '+str(row[4]))
                                job.meta['err_msg'] = 'rr_valueerror'
                                job.save_meta()
                                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                                #return False, 'rr_valueerror', b''
                                return False, 'rr_valueerror', ''
                shutil.copy2(outfile_path, rules_file)
            ####################################
            ########## Retropath2 ##############
            ####################################
            logger.info('---------- RP2 -----')
            job.meta['progress'] = 5
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            rp2_file = os.path.join(rp_results_folder, 'rp2_pathways.csv')
            logger.debug('Rules file: '+str(rules_file))
            logger.debug('Timeout: '+str(sub_timeout*60.0)+' seconds')
            dmin = 0
            dmax = 1000
            mwmax_source = 1000
            mwmax_cof = 1000
            logger.debug('Default RP2 configrations:')
            logger.debug('dmin: '+str(dmin))
            logger.debug('dmanx: '+str(dmax))
            logger.debug('mwmax_source: '+str(mwmax_source))
            logger.debug('mwmax_cof: '+str(mwmax_cof))
            is_timeout = False
            is_results_empty = True
            ### run the KNIME RETROPATH2.0 workflow
            with tempfile.TemporaryDirectory() as tmp_rp2_folder:
                knime_command = KPATH+' -nosplash -nosave -reset --launcher.suppressErrors -application org.knime.product.KNIME_BATCH_APPLICATION -workflowFile='+RP_WORK_PATH+' -workflow.variable=input.dmin,"'+str(dmin)+'",int -workflow.variable=input.dmax,"'+str(dmax)+'",int -workflow.variable=input.max-steps,"'+str(max_steps)+'",int -workflow.variable=input.sourcefile,"'+str(source_file)+'",String -workflow.variable=input.sinkfile,"'+str(sink_file)+'",String -workflow.variable=input.rulesfile,"'+str(rules_file)+'",String -workflow.variable=input.topx,"'+str(topx)+'",int -workflow.variable=input.mwmax-source,"'+str(mwmax_source)+'",int -workflow.variable=input.mwmax-cof,"'+str(mwmax_cof)+'",int -workflow.variable=output.dir,"'+str(tmp_rp2_folder)+'/",String -workflow.variable=output.solutionfile,"results.csv",String -workflow.variable=output.sourceinsinkfile,"source-in-sink.csv",String'
                logger.debug('KNIME command: '+str(knime_command))
                commandObj = subprocess.Popen(knime_command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=limit_virtual_memory)
                try:
                    stdout = b''
                    stderr = b''
                    try:
                        logger.debug('Running...')
                        stdout, stderr = commandObj.communicate(timeout=sub_timeout*60.0) #subprocess timeout is in seconds while we input minutes
                        logger.debug('Ran RetroPath2.0!')
                        stdout = stdout.decode('utf-8')
                        stderr = stderr.decode('utf-8')
                    except subprocess.TimeoutExpired as e:
                        logger.warning('RetroPath2.0 has reached its execution timeout limit')
                        commandObj.kill()
                        is_timeout = True
                    except subprocess.CalledProcessError as e:
                        logger.error('Non-zero raised by RetroPath2.0')
                        logger.error(e)
                        logger.debug(stdout)
                        logger.debug(stderr)
                    except subprocess.SubprocessError as e:
                        logger.error('The subprocess throws an error')
                        logger.error(e)
                    except:
                        logger.error('All other exceptions')
                        e = sys.exc_info()[0]
                        logger.error(e)
                    logger.debug('Output folder: '+str(glob.glob(os.path.join(tmp_rp2_folder, '*'))))
                    #check to see if the results.csv is empty
                    logger.debug('Checking the results.csv file')
                    try:
                        count = 0
                        with open(os.path.join(tmp_rp2_folder, 'results.csv')) as f:
                            reader = csv.reader(f, delimiter=',', quotechar='"')
                            for i in reader:
                                count += 1
                        if count>1:
                            is_results_empty = False
                    except (IndexError, FileNotFoundError) as e:
                        logger.debug('No results.csv file')
                        pass
                    ##################### HANDLE all the different cases ###################
                    ### if source is in sink. Note making sure that it contains more than the default first line
                    ### if java has an memory issue
                    if 'There is insufficient memory for the Java Runtime Environment to continue' in stdout:
                        if not is_results_empty and partial_retro:
                            logger.warning('RetroPath2.0 does not have sufficient memory to continue')
                            shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                            logger.warning('Passing the results file instead')
                        else:
                            logger.error('RetroPath2.0 does not have sufficient memory to continue')
                            job.meta['err_msg'] = 'rp2_ram'
                            job.save_meta()
                            modResJSON(job.id, job_meta=job.meta, job_status='failed')
                            #return False, 'rp2_ram', b''
                            return False, 'rp2_ram', ''
                    ### handle timeout
                    if is_timeout:
                        if not is_results_empty and partial_retro:
                            logger.warning('Timeout from retropath2.0 ('+str(sub_timeout)+' minutes)')
                            shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                        else:
                            logger.error('Timeout from retropath2.0 ('+str(sub_timeout)+' minutes)')
                            job.meta['err_msg'] = 'rp2_timeout'
                            job.save_meta()
                            modResJSON(job.id, job_meta=job.meta, job_status='failed')
                            #return False, 'rp2_timeout', b''
                            return False, 'rp2_timeout', ''
                    try:
                        count = 0
                        with open(os.path.join(tmp_rp2_folder, 'source-in-sink.csv')) as f:
                            reader = csv.reader(f, delimiter=',', quotechar='"')
                            for i in reader:
                                count += 1
                        if count>1:
                            logger.error('Execution problem of RetroPath2.0. Source has been found in the sink')
                            job.meta['err_msg'] = 'rp2_sourceinsink'
                            job.save_meta()
                            modResJSON(job.id, job_meta=job.meta, job_status='failed')
                            #return False, 'rp2_sourceinsink', b''
                            return False, 'rp2_sourceinsink', ''
                    except FileNotFoundError as e:
                        logger.error('Cannot find source-in-sink.csv file. Probably an execution error.')
                        logger.error(e)
                        job.meta['err_msg'] = 'rp2_execerror'
                        job.save_meta()
                        modResJSON(job.id, job_meta=job.meta, job_status='failed')
                        #return False, 'rp2_execerror', b''
                        return False, 'rp2_execerror', ''
                    ############## IF ALL IS GOOD ##############
                    ### csv scope copy to the .dat location
                    try:
                        csv_scope = glob.glob(os.path.join(tmp_rp2_folder, '*_scope.csv'))
                        shutil.copy2(csv_scope[0], rp2_file)
                    except IndexError as e:
                        if not is_results_empty and partial_retro:
                            logger.warning('No scope file generated')
                            shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                            logger.warning('Passing the results file instead')
                        else:
                            logger.error('RetroPath2.0 has not found any results')
                            job.meta['err_msg'] = 'rp2_noresults'
                            job.save_meta()
                            modResJSON(job.id, job_meta=job.meta, job_status='failed')
                            #return False, 'rp2_noresults', b''
                            return False, 'rp2_noresults', ''
                except OSError as e:
                    if not is_results_empty and partial_retro:
                        logger.warning('Running the RetroPath2.0 Knime program produced an OSError')
                        logger.warning(e)
                        shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                        logger.warning('Passing the results file instead')
                    else:
                        logger.error('Running the RetroPath2.0 Knime program produced an OSError')
                        logger.error(e)
                        job.meta['err_msg'] = 'rp2_oserror'
                        job.save_meta()
                        modResJSON(job.id, job_meta=job.meta, job_status='failed')
                        #return False, 'rp2_oserror', b''
                        return False, 'rp2_oserror', ''
                except ValueError as e:
                    if not is_results_empty and partial_retro:
                        logger.warning('Cannot set the RAM usage limit')
                        logger.warning(e)
                        shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                        logger.warning('Passing the results file instead')
                    else:
                        logger.error('Cannot set the RAM usage limit')
                        logger.error(e)
                        job.meta['err_msg'] = 'rp2_ram'
                        job.save_meta()
                        modResJSON(job.id, job_meta=job.meta, job_status='failed')
                        #return False, 'rp2_ram', b''
                        return False, 'rp2_ram', ''
            ###################################
            ######### rp2paths ################
            ###################################
            logger.info('---------- RP2paths -----')
            job.meta['progress'] = 6
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            rp2paths_pathways_file = os.path.join(rp_results_folder, 'rp2paths_pathways.csv')
            rp2paths_compounds_file = os.path.join(rp_results_folder, 'rp2paths_compounds.tsv')
            with tempfile.TemporaryDirectory() as tmp_rp2paths_folder:
                rp2paths_command = 'python /home/rp2paths/RP2paths.py all '+str(rp2_file)+' --outdir '+str(tmp_rp2paths_folder)+' --timeout '+str(int(sub_timeout*60.0))
                try:
                    commandObj = subprocess.Popen(rp2paths_command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=limit_virtual_memory)
                    result = b''
                    error = b''
                    result, error = commandObj.communicate()
                    result = result.decode('utf-8')
                    error = error.decode('utf-8')
                    #TODO test to see what is the correct phrase
                    if 'TIMEOUT' in result:
                        logger.error('Timeout from of ('+str(timeout)+' minutes)')
                        job.meta['err_msg'] = 'rp2paths_timeout'
                        job.save_meta()
                        modResJSON(job.id, job_meta=job.meta, job_status='failed')
                        #return False, 'rp2paths_timeout', b''
                        return False, 'rp2paths_timeout', ''
                    if 'failed to map segment from shared object' in error:
                        logger.error('RP2paths does not have sufficient memory to continue')
                        job.meta['err_msg'] = 'rp2paths_ram'
                        job.save_meta()
                        modResJSON(job.id, job_meta=job.meta, job_status='failed')
                        #return False, 'rp2paths_ram', b''
                        return False, 'rp2paths_ram', ''
                    ### convert the result to binary and return ###
                    logger.debug(glob.glob(os.path.join(tmp_rp2paths_folder, '*')))
                    try:
                        shutil.copy2(os.path.join(tmp_rp2paths_folder, 'out_paths.csv'), rp2paths_pathways_file)
                        shutil.copy2(os.path.join(tmp_rp2paths_folder, 'compounds.txt'), rp2paths_compounds_file)
                    except FileNotFoundError as e:
                        logger.error('Cannot find the output files out_paths.csv or compounds.txt')
                        job.meta['err_msg'] = 'rp2paths_outputerr'
                        job.save_meta()
                        modResJSON(job.id, job_meta=job.meta, job_status='failed')
                        #return False, 'rp2paths_outputerr', b''
                        return False, 'rp2paths_outputerr', ''
                except OSError as e:
                    logger.error('Subprocess detected an error when calling the rp2paths command')
                    job.meta['err_msg'] = 'rp2paths_subprocess'
                    job.save_meta()
                    modResJSON(job.id, job_meta=job.meta, job_status='failed')
                    #return False, 'rp2paths_subprocess', b''
                    return False, 'rp2paths_subprocess', ''
                except ValueError as e:
                    logger.error('Cannot set the RAM usage limit')
                    job.meta['err_msg'] = 'rp2paths_ram'
                    job.save_meta()
                    modResJSON(job.id, job_meta=job.meta, job_status='failed')
                    #return False, 'rp2paths_ram', b''
                    return False, 'rp2paths_ram', ''
            ##################################
            ####### Analysis #################
            ##################################
            logger.info('---------- rpReader -----')
            job.meta['progress'] = 7
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            rpreader_status = rpReader.rp2ToCollection(rp2_file,
                                                       rp2paths_compounds_file,
                                                       rp2paths_pathways_file,
                                                       rpcollection_file,
                                                       rpcache=global_rpcache)
            if not rpreader_status:
                logger.error('Problem running rpReader')
                job.meta['err_msg'] = 'rpreader'
                job.save_meta()
                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                #return False, 'rpreader', b''
                return False, 'rpreader', ''
            logger.info('---------- rpEquilibrator -----')
            job.meta['progress'] = 8
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            logger.debug('running Eq')
            rpeq_status = rpEquilibrator.runCollection(rpcollection_file,
                                                       rpcollection_file,
                                                       ph=ph,
                                                       ionic_strength=ionic_strength,
                                                       temp_k=temp_k,
                                                       rpcache=global_rpcache)
            if not rpeq_status:
                logger.error('Problem running rpEquilibrator')
                job.meta['err_msg'] = 'rpeq'
                job.save_meta()
                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                #return False, 'rpeq', b''
                return False, 'rpeq', ''
            logger.info('---------- rpFBA -----')
            job.meta['progress'] = 9
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            logger.debug('running FBA')
            rpfba_status = rpFBA.runCollection(rpcollection_file,
                                               gem_file,
                                               rpcollection_file,
                                               num_workers=1,
                                               keep_merged=True,
                                               del_sp_pro=False,
                                               del_sp_react=False,
                                               rpcache=global_rpcache)
            if not rpfba_status:
                logger.error('Problem running rpFBA')
                job.meta['err_msg'] = 'rpfba'
                job.save_meta()
                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                #return False, 'rpfba', b''
                return False, 'rpfba', ''
            logger.info('---------- rpSelenzyme -----')
            job.meta['progress'] = 10
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            rpsel_status = rpSelenzyme.runCollection(rpcollection_file,
                                                     taxo_id,
                                                     rpcollection_file,
                                                     cache_path=os.path.join('/home', 'metaxime', 'input_cache', 'rpselenzyme_data.tar.xz'),
                                                     rpcache=global_rpcache)
            if not rpsel_status:
                logger.error('Problem running rpSelenzyme')
                job.meta['err_msg'] = 'rpsel'
                job.save_meta()
                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                #return False, 'rpsel', b''
                return False, 'rpsel', ''
            logger.info('---------- rpGlobalScore -----')
            job.meta['progress'] = 11
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            rpglo_status = rpGlobalScore.runCollection(rpcollection_file,
                                                       rpcollection_file,
                                                       rpcache=global_rpcache)
            if not rpglo_status:
                logger.error('Problem running rpGlobalScore')
                job.meta['err_msg'] = 'rpglo'
                job.save_meta()
                modResJSON(job.id, job_meta=job.meta, job_status='failed')
                #return False, 'rpglo', b''
                return False, 'rpglo', ''
            job.meta['progress'] = 12
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            ################ generate the various netowrk JSON and SBML JSON #####
            logger.info('------------ JSON ----------')
            rpSBML.batchAsDict(rpcollection_file)
            rpGraph.batchNetworkDict(rpcollection_file, is_gem_sbml=False)
            ################ make json with the rank and the filenames #########
            logger.info('----------- save results -------')
            #copy the results to the shared folder
            result_output_path = os.path.join('/mx-results/', str(job.id))
            shutil.move(tmp_output_folder, result_output_path)
            subprocess.call(['tar', '-xf', os.path.join(result_output_path, 'rpcollection.tar.xz'), '-C', result_output_path])
            subprocess.call(['chmod', '-R', '777', result_output_path])
            job.meta['progress'] = 13
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta)
            #### update the sttus json ###
            modResJSON(job.id,
                       job_status='finished',
                       job_meta=job.meta,
                       output_path=result_output_path,
                       output_tar=os.path.join(result_output_path, 'rpcollection.tar.xz'))
            logger.info('-------- done --------')
            return True, 'success', rpcollection_file
        except:
            logger.error("Unexpected error:", sys.exc_info()[0])
            job.meta['err_msg'] = 'unexpected'
            job.save_meta()
            modResJSON(job.id, job_meta=job.meta, job_status='failed')
            return False, 'failure', sys.exc_info()[0]
