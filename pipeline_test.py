##############################################################
########################## PIPELINE ##########################
##############################################################

import tempfile
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
from equilibrator_api import ComponentContribution


from metaxime import rpCache
from metaxime import rpGraph
from metaxime import rpReader
from metaxime import rpFBA
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

print(__name__)


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
    print('####### pipieline ######')
    #modResJSON(job.id, job_meta=job.meta, job_status='running')
    print('cache')
    global_rpcache = rpCache()
    print('populate cache')
    global_rpcache.populateCache()
    print(target_smiles)
    print(gem_name)
    print(max_steps)
    target_inchi = MolToInchi(MolFromSmiles(target_smiles, sanitize=True))
    print(target_inchi)
    with tempfile.TemporaryDirectory() as tmp_output_folder:
        ############# source file #############
        print('------ source file -----')
        rpcollection_file = os.path.join(tmp_output_folder, 'rpcollection.tar.xz')
        rp_results_folder = os.path.join(tmp_output_folder, 'rp_results')
        os.mkdir(rp_results_folder)
        source_file = os.path.join(rp_results_folder, 'source.csv')
        with open(source_file, 'w') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(['Name', 'InChI'])
            filewriter.writerow(['target', target_inchi])
        ############ sink file ##############
        print('-------- sink file --------')
        try:
            gem_file = model_list[gem_name]
            taxo_id = model_taxo_id[gem_name]
            sink_file = sink_list[gem_name]
        except KeyError:
            print(gem_name)
            #return False, 'gem', b''
            return False, 'gem', ''
        '''#only when passing an SBML file to process
        sink_file = os.path.join(tmp_output_folder, 'sink.csv')
        gensink_status = rpExtractSink.genSink(gem_file, sink_file, remove_dead_end=True)
        if not gensink_status:
            print('Problem generating sink')
            return False, 'gensink'
        '''
        shutil.copy2(gem_file, os.path.join(rp_results_folder, 'gem_sbml.xml'))
        shutil.copy2(sink_file, os.path.join(rp_results_folder, 'sink.csv'))
        ####################################
        ########## Reaction Rules ##########
        ####################################
        print('---------- reaction rules -----')
        rules_file = os.path.join(rp_results_folder, 'reaction_rules.csv')
        rule_file = None
        if rules_type=='all':
            rule_file = '/home/retrorules/rules_rall_rp2.csv'
        elif rules_type=='forward':
            rule_file = '/home/retrorules/rules_rall_rp2_forward.csv'
        elif rules_type=='retro':
            rule_file = '/home/retrorules/rules_rall_rp2_retro.csv'
        else:
            print(rules_type)
            #return False, 'rr_inputerror', b''
            return False, 'rr_inputerror', ''
        #check the input diameters are valid #
        try:
            s_diameters = [int(i) for i in rules_diameters.split(',')]
            valid_diameters = []
            for i in s_diameters:
                if i not in [2,4,6,8,10,12,14,16]:
                    print(i)
                else:
                    valid_diameters.append(i)
        except ValueError:
            print('RR: Invalid diamter entry. Must be int of either 2,4,6,8,10,12,14,16')
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
                            print(row[4])
                            #return False, 'rr_valueerror', b''
                            return False, 'rr_valueerror', ''
            shutil.copy2(outfile_path, rules_file)
        ####################################
        ########## Retropath2 ##############
        ####################################
        print('---------- RP2 -----')
        rp2_file = os.path.join(rp_results_folder, 'rp2_pathways.csv')
        print(rules_file)
        dmin = 0
        dmax = 1000
        mwmax_source = 1000
        mwmax_cof = 1000
        print('Default RP2 configrations:')
        print(dmin)
        print(dmax)
        print(mwmax_source)
        print(mwmax_cof)
        is_timeout = False
        is_results_empty = True
        ### run the KNIME RETROPATH2.0 workflow
        with tempfile.TemporaryDirectory() as tmp_rp2_folder:
            knime_command = KPATH+' -nosplash -nosave -reset --launcher.suppressErrors -application org.knime.product.KNIME_BATCH_APPLICATION -workflowFile='+RP_WORK_PATH+' -workflow.variable=input.dmin,"'+str(dmin)+'",int -workflow.variable=input.dmax,"'+str(dmax)+'",int -workflow.variable=input.max-steps,"'+str(max_steps)+'",int -workflow.variable=input.sourcefile,"'+str(source_file)+'",String -workflow.variable=input.sinkfile,"'+str(sink_file)+'",String -workflow.variable=input.rulesfile,"'+str(rules_file)+'",String -workflow.variable=input.topx,"'+str(topx)+'",int -workflow.variable=input.mwmax-source,"'+str(mwmax_source)+'",int -workflow.variable=input.mwmax-cof,"'+str(mwmax_cof)+'",int -workflow.variable=output.dir,"'+str(tmp_rp2_folder)+'/",String -workflow.variable=output.solutionfile,"results.csv",String -workflow.variable=output.sourceinsinkfile,"source-in-sink.csv",String'
            print(knime_command)
            commandObj = subprocess.Popen(knime_command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=limit_virtual_memory)
            try:
                stdout = b''
                stderr = b''
                try:
                    print('Running...')
                    stdout, stderr = commandObj.communicate(timeout=sub_timeout*60.0) #subprocess timeout is in seconds while we input minutes
                    print('Ran RetroPath2.0!')
                    stdout = stdout.decode('utf-8')
                    stderr = stderr.decode('utf-8')
                except subprocess.TimeoutExpired as e:
                    print('RetroPath2.0 has reached its execution timeout limit')
                    commandObj.kill()
                    is_timeout = True
                except subprocess.CalledProcessError as e:
                    print('Non-zero raised by RetroPath2.0')
                    print(e)
                    print(stdout)
                    print(stderr)
                except subprocess.SubprocessError as e:
                    print('The subprocess throws an error')
                    print(e)
                except:
                    print('All other exceptions')
                    e = sys.exc_info()[0]
                    print(e)
                print(tmp_rp2_folder, '*')
                #check to see if the results.csv is empty
                print('Checking the results.csv file')
                try:
                    count = 0
                    with open(os.path.join(tmp_rp2_folder, 'results.csv')) as f:
                        reader = csv.reader(f, delimiter=',', quotechar='"')
                        for i in reader:
                            count += 1
                    if count>1:
                        is_results_empty = False
                except (IndexError, FileNotFoundError) as e:
                    print('No results.csv file')
                    pass
                ##################### HANDLE all the different cases ###################
                ### if source is in sink. Note making sure that it contains more than the default first line
                ### if java has an memory issue
                if 'There is insufficient memory for the Java Runtime Environment to continue' in stdout:
                    if not is_results_empty and partial_retro:
                        print('RetroPath2.0 does not have sufficient memory to continue')
                        shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                        print('Passing the results file instead')
                    else:
                        print('RetroPath2.0 does not have sufficient memory to continue')
                        #return False, 'rp2_ram', b''
                        return False, 'rp2_ram', ''
                ### handle timeout
                if is_timeout:
                    if not is_results_empty and partial_retro:
                        shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                    else:
                        #return False, 'rp2_timeout', b''
                        return False, 'rp2_timeout', ''
                try:
                    count = 0
                    with open(os.path.join(tmp_rp2_folder, 'source-in-sink.csv')) as f:
                        reader = csv.reader(f, delimiter=',', quotechar='"')
                        for i in reader:
                            count += 1
                    if count>1:
                        print('Execution problem of RetroPath2.0. Source has been found in the sink')
                        #return False, 'rp2_sourceinsink', b''
                        return False, 'rp2_sourceinsink', ''
                except FileNotFoundError as e:
                    print('Cannot find source-in-sink.csv file. Probably an execution error.')
                    print(e)
                    #return False, 'rp2_execerror', b''
                    return False, 'rp2_execerror', ''
                ############## IF ALL IS GOOD ##############
                ### csv scope copy to the .dat location
                try:
                    csv_scope = glob.glob(os.path.join(tmp_rp2_folder, '*_scope.csv'))
                    shutil.copy2(csv_scope[0], rp2_file)
                except IndexError as e:
                    if not is_results_empty and partial_retro:
                        print('No scope file generated')
                        shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                        print('Passing the results file instead')
                    else:
                        print('RetroPath2.0 has not found any results')
                        #return False, 'rp2_noresults', b''
                        return False, 'rp2_noresults', ''
            except OSError as e:
                if not is_results_empty and partial_retro:
                    print('Running the RetroPath2.0 Knime program produced an OSError')
                    print(e)
                    shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                    print('Passing the results file instead')
                else:
                    print('Running the RetroPath2.0 Knime program produced an OSError')
                    print(e)
                    #return False, 'rp2_oserror', b''
                    return False, 'rp2_oserror', ''
            except ValueError as e:
                if not is_results_empty and partial_retro:
                    print('Cannot set the RAM usage limit')
                    print(e)
                    shutil.copy2(os.path.join(tmp_rp2_folder, 'results.csv'), rp2_file)
                    print('Passing the results file instead')
                else:
                    print('Cannot set the RAM usage limit')
                    print(e)
                    #return False, 'rp2_ram', b''
                    return False, 'rp2_ram', ''
        ###################################
        ######### rp2paths ################
        ###################################
        print('---------- RP2paths -----')
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
                    #return False, 'rp2paths_timeout', b''
                    return False, 'rp2paths_timeout', ''
                if 'failed to map segment from shared object' in error:
                    print('RP2paths does not have sufficient memory to continue')
                    #return False, 'rp2paths_ram', b''
                    return False, 'rp2paths_ram', ''
                ### convert the result to binary and return ###
                print(tmp_rp2paths_folder, '*')
                try:
                    shutil.copy2(os.path.join(tmp_rp2paths_folder, 'out_paths.csv'), rp2paths_pathways_file)
                    shutil.copy2(os.path.join(tmp_rp2paths_folder, 'compounds.txt'), rp2paths_compounds_file)
                except FileNotFoundError as e:
                    print('Cannot find the output files out_paths.csv or compounds.txt')
                    #return False, 'rp2paths_outputerr', b''
                    return False, 'rp2paths_outputerr', ''
            except OSError as e:
                print('Subprocess detected an error when calling the rp2paths command')
                #return False, 'rp2paths_subprocess', b''
                return False, 'rp2paths_subprocess', ''
            except ValueError as e:
                print('Cannot set the RAM usage limit')
                #return False, 'rp2paths_ram', b''
                return False, 'rp2paths_ram', ''
        ##################################
        ####### Analysis #################
        ##################################
        print('---------- rpReader -----')
        print(glob.glob(os.path.join(tmp_output_folder, '*')))
        print(rpcollection_file)
        rpreader_status = rpReader.rp2ToCollection(rp2_file,
                                                   rp2paths_compounds_file,
                                                   rp2paths_pathways_file,
                                                   rpcollection_file,
                                                   rpcache=global_rpcache)
        if not rpreader_status:
            print('Problem running rpReader')
            #return False, 'rpreader', b''
            return False, 'rpreader', ''
        print('---------- rpEquilibrator -----')
        print(rpcollection_file)
        print(glob.glob(os.path.join(tmp_output_folder, '*')))
        print('running Eq')
        rpeq_status = rpEquilibrator.runCollection(rpcollection_file,
                                                   rpcollection_file,
                                                   ph=ph,
                                                   ionic_strength=ionic_strength,
                                                   temp_k=temp_k,
                                                   rpcache=global_rpcache)
        if not rpeq_status:
            print('Problem running rpEquilibrator')
            #return False, 'rpeq', b''
            return False, 'rpeq', ''
        print('---------- rpFBA -----')
        print(glob.glob(os.path.join(tmp_output_folder, '*')))
        print('running FBA')
        rpfba_status = rpFBA.runCollection(rpcollection_file,
                                           gem_file,
                                           rpcollection_file,
                                           num_workers=1,
                                           keep_merged=True,
                                           del_sp_pro=False,
                                           del_sp_react=False,
                                           rpcache=global_rpcache)
        if not rpfba_status:
            print('Problem running rpFBA')
            #return False, 'rpfba', b''
            return False, 'rpfba', ''
        print('---------- rpSelenzyme -----')
        rpsel_status = rpSelenzyme.runCollection(rpcollection_file,
                                                 taxo_id,
                                                 rpcollection_file,
                                                 cache_path=os.path.join('/home', 'metaxime', 'input_cache', 'rpselenzyme_data.tar.xz'),
                                                 rpcache=global_rpcache)
        if not rpsel_status:
            print('Problem running rpSelenzyme')
            #return False, 'rpsel', b''
            return False, 'rpsel', ''
        print('---------- rpGlobalScore -----')
        rpglo_status = rpGlobalScore.runCollection(rpcollection_file,
                                                   rpcollection_file,
                                                   rpcache=global_rpcache)
        if not rpglo_status:
            print('Problem running rpGlobalScore')
            #return False, 'rpglo', b''
            return False, 'rpglo', ''
        ################ generate the various netowrk JSON and SBML JSON #####
        print('------------ JSON ----------')
        rpSBML.batchAsDict(rpcollection_file)
        rpGraph.batchNetworkDict(rpcollection_file, is_gem_sbml=False)
        ################ make json with the rank and the filenames #########
        print('----------- save results -------')
        #copy the results to the shared folder
        result_output_path = os.path.join('/mx-results/', 'test')
        shutil.move(tmp_output_folder, result_output_path)
        subprocess.call(['chmod', '-R', '777', result_output_path])
        #### update the sttus json ###
        print('-------- done --------')
        return True, 'success', rpcollection_file
