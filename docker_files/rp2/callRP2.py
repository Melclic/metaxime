#!/usr/bin/env python3
"""
Created on January 16 2020

@author: Melchior du Lac
@description: Galaxy script to query rpRetroPath2.0 REST service

"""
import subprocess
import logging
import csv
import glob
import resource
import tempfile
import sys
import argparse
import tarfile
import shutil
import os


KPATH = '/usr/local/knime/knime'
RP_WORK_PATH = '/home/RetroPath2.0.knwf'


logging.basicConfig(
    #level=logging.DEBUG,
    #level=logging.WARNING,
    level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


#MAX_VIRTUAL_MEMORY = 20000*1024*1024 # 20 GB -- define what is the best
MAX_VIRTUAL_MEMORY = 30000*1024*1024 # 30 GB -- define what is the best

def limit_virtual_memory():
    """Limit the virtual of the subprocess call
    """
    resource.setrlimit(resource.RLIMIT_AS, (MAX_VIRTUAL_MEMORY, resource.RLIM_INFINITY))


def run_rp2(source_path, sink_path, rules_path, max_steps, topx=100, dmin=0, dmax=1000, mwmax_source=1000, mwmax_cof=1000, timeout=30, partial_retro=False, logger=None):
    """Call the KNIME RetroPath2.0 workflow

    :param source_path: The path to the source file
    :param sink_path: The path to the sink file
    :param rules_path: The path to the rules file
    :param max_steps: The maximal number of steps
    :param topx: The top number of reaction rules to keep at each iteraction (Default: 100)
    :param dmin: The minimum diameter of the reaction rules (Default: 0)
    :param dmax: The miximum diameter of the reaction rules (Default: 1000)
    :param mwmax_source: The maximal molecular weight of the intermediate compound (Default: 1000)
    :param mwmax_cof: The coefficient of the molecular weight of the intermediate compound (Default: 1000)
    :param timeout: The timeout of the function in minutes (Default: 30)
    :param ram_limit: Set the upper bound of the RAM usage of the tool (Default: 30 GB)
    :param partial_retro: Return partial results if the execution is interrupted for any reason (Default: False)
    :param logger: Logger object (Default: None)

    :type source_path: str
    :type sink_path: str
    :type rules_path: str
    :type max_steps: int
    :type topx: int
    :type dmin: int
    :type dmax: int
    :type mwmax_source: int
    :type mwmax_cof: int
    :type timeout: int
    :type partial_retro: bool
    :type logger: logging

    :rtype: tuple
    :return: tuple of bytes with the results, the status message, the KNIME command used
    """
    if logger==None:
        logger = logging.getLogger(__name__)
    logger.debug('Rules file: '+str(rules_path))
    logger.debug('Timeout: '+str(timeout*60.0)+' seconds')
    is_timeout = False
    is_results_empty = True
    ### run the KNIME RETROPATH2.0 workflow
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        try:
            knime_command = KPATH+' -nosplash -nosave -reset --launcher.suppressErrors -application org.knime.product.KNIME_BATCH_APPLICATION -workflowFile='+RP_WORK_PATH+' -workflow.variable=input.dmin,"'+str(dmin)+'",int -workflow.variable=input.dmax,"'+str(dmax)+'",int -workflow.variable=input.max-steps,"'+str(max_steps)+'",int -workflow.variable=input.sourcefile,"'+str(source_path)+'",String -workflow.variable=input.sinkfile,"'+str(sink_path)+'",String -workflow.variable=input.rulesfile,"'+str(rules_path)+'",String -workflow.variable=input.topx,"'+str(topx)+'",int -workflow.variable=input.mwmax-source,"'+str(mwmax_source)+'",int -workflow.variable=input.mwmax-cof,"'+str(mwmax_cof)+'",int -workflow.variable=output.dir,"'+str(tmpOutputFolder)+'/",String -workflow.variable=output.solutionfile,"results.csv",String -workflow.variable=output.sourceinsinkfile,"source-in-sink.csv",String'
            logger.debug('KNIME command: '+str(knime_command))
            commandObj = subprocess.Popen(knime_command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=limit_virtual_memory)
            result = b''
            error = b''
            try:
                #commandObj.wait(timeout=timeout) #subprocess timeout is in seconds while we input minutes
                result, error = commandObj.communicate(timeout=timeout*60.0) #subprocess timeout is in seconds while we input minutes
            except subprocess.TimeoutExpired as e:
                logging.warning('RetroPath2.0 has reached its execution timeout limit')
                commandObj.kill()
                is_timeout = True
            #(result, error) = commandObj.communicate()
            result = result.decode('utf-8')
            error = error.decode('utf-8')
            logger.debug('RetroPath2.0 results message: '+str(result))
            logger.debug('RetroPath2.0 error message: '+str(error))
            logger.debug('Output folder: '+str(glob.glob(tmpOutputFolder+'/*')))
            #check to see if the results.csv is empty
            try:
                count = 0
                with open(tmpOutputFolder+'/results.csv') as f:
                    reader = csv.reader(f, delimiter=',', quotechar='"')
                    for i in reader:
                        count += 1
                if count>1:
                    is_results_empty = False
            except (IndexError, FileNotFoundError) as e:
                logger.debug('No results.csv file')
                pass
            ########################################################################
            ##################### HANDLE all the different cases ###################
            ########################################################################
            ### if source is in sink. Note making sure that it contains more than the default first line
            try:
                count = 0
                with open(tmpOutputFolder+'/source-in-sink.csv') as f:
                    reader = csv.reader(f, delimiter=',', quotechar='"')
                    for i in reader:
                        count += 1
                if count>1:
                    logger.error('Source has been found in the sink')
                    return b'', b'sourceinsinkerror', str('Command: '+str(knime_command)+'\n Error: Source found in sink\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
            except FileNotFoundError as e:
                logger.error('Cannot find source-in-sink.csv file')
                logger.error(e)
                return b'', b'sourceinsinknotfounderror', str('Command: '+str(knime_command)+'\n Error: '+str(e)+'\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
            ### handle timeout
            if is_timeout:
                if not is_results_empty and partial_retro:
                    logger.warning('Timeout from retropath2.0 ('+str(timeout)+' minutes)')
                    with open(tmpOutputFolder+'/results.csv', 'rb') as op:
                        results_csv = op.read()
                    return results_csv, b'timeoutwarning', str('Command: '+str(knime_command)+'\n Error: Memory error \n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
                else:
                    logger.error('Timeout from retropath2.0 ('+str(timeout)+' minutes)')
                    return b'', b'timeouterror', str('Command: '+str(knime_command)+'\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
            ### if java has an memory issue
            if 'There is insufficient memory for the Java Runtime Environment to continue' in result:
                if not is_results_empty and partial_retro:
                    logger.warning('RetroPath2.0 does not have sufficient memory to continue')
                    with open(tmpOutputFolder+'/results.csv', 'rb') as op:
                        results_csv = op.read()
                    logger.warning('Passing the results file instead')
                    return results_csv, b'memwarning', str('Command: '+str(knime_command)+'\n Error: Memory error \n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
                else:
                    logger.error('RetroPath2.0 does not have sufficient memory to continue')
                    return b'', b'memerror', str('Command: '+str(knime_command)+'\n Error: Memory error \n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
            ############## IF ALL IS GOOD ##############
            ### csv scope copy to the .dat location
            try:
                csv_scope = glob.glob(tmpOutputFolder+'/*_scope.csv')
                with open(csv_scope[0], 'rb') as op:
                    scope_csv = op.read()
                return scope_csv, b'noerror', str('').encode('utf-8')
            except IndexError as e:
                if not is_results_empty and partial_retro:
                    logger.warning('No scope file generated')
                    with open(tmpOutputFolder+'/results.csv', 'rb') as op:
                        results_csv = op.read()
                    logger.warning('Passing the results file instead')
                    return results_csv, b'noresultwarning', str('Command: '+str(knime_command)+'\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
                else:
                    logger.error('RetroPath2.0 has not found any results')
                    return b'', b'noresulterror', str('Command: '+str(knime_command)+'\n Error: '+str(e)+'\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
        except OSError as e:
            if not is_results_empty and partial_retro:
                logger.warning('Running the RetroPath2.0 Knime program produced an OSError')
                logger.warning(e) 
                with open(tmpOutputFolder+'/results.csv', 'rb') as op:
                    results_csv = op.read()
                logger.warning('Passing the results file instead')
                return results_csv, b'oswarning', str('Command: '+str(knime_command)+'\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
            else:
                logger.error('Running the RetroPath2.0 Knime program produced an OSError')
                logger.error(e)
                return b'', b'oserror', str('Command: '+str(knime_command)+'\n Error: '+str(e)+'\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
        except ValueError as e:
            if not is_results_empty and partial_retro:
                logger.warning('Cannot set the RAM usage limit')
                logger.warning(e)
                with open(tmpOutputFolder+'/results.csv', 'rb') as op:
                    results_csv = op.read()
                logger.warning('Passing the results file instead')
                return results_csv, b'ramwarning', str('Command: '+str(knime_command)+'\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')
            else:
                logger.error('Cannot set the RAM usage limit')
                logger.error(e)
                return b'', b'ramerror', str('Command: '+str(knime_command)+'\n Error: '+str(e)+'\n tmpOutputFolder: '+str(glob.glob(tmpOutputFolder+'/*'))).encode('utf-8')


if __name__ == "__main__":
    #### WARNING: as it stands one can only have a single source molecule
    parser = argparse.ArgumentParser('Python wrapper for the KNIME workflow to run RetroPath2.0')
    parser.add_argument('-sinkfile', type=str)
    parser.add_argument('-sourcefile', type=str)
    parser.add_argument('-max_steps', type=int)
    parser.add_argument('-rulesfile', type=str)
    parser.add_argument('-rulesfile_format', type=str)
    parser.add_argument('-scope_csv', type=str)
    parser.add_argument('-topx', type=int, default=100)
    parser.add_argument('-dmin', type=int, default=0)
    parser.add_argument('-dmax', type=int, default=100)
    parser.add_argument('-mwmax_source', type=int, default=1000)
    parser.add_argument('-mwmax_cof', type=int, default=1000)
    parser.add_argument('-timeout', type=int, default=90)
    parser.add_argument('-partial_retro', type=str, default='False')
    params = parser.parse_args()
    if params.max_steps<=0:
        logging.error('Maximal number of steps cannot be less or equal to 0')
        exit(1)
    if params.topx<0:
        logging.error('Cannot have a topx value that is <0: '+str(params.topx))
        exit(1)
    if params.dmin<0:
        logging.error('Cannot have a dmin value that is <0: '+str(params.dmin))
        exit(1)
    if params.dmax<0:
        logging.error('Cannot have a dmax value that is <0: '+str(params.dmax))
        exit(1)
    if params.dmax>1000:
        logging.error('Cannot have a dmax valie that is >1000: '+str(params.dmax))
        exit(1)
    if params.dmax<params.dmin:
        logging.error('Cannot have dmin>dmax : dmin: '+str(params.dmin)+', dmax: '+str(params.dmax))
        exit(1)
    if params.partial_retro=='False' or params.partial_retro=='false' or params.partial_retro=='F':
        partial_retro = False
    elif params.partial_retro=='True' or params.partial_retro=='true' or params.partial_retro=='T':
        partial_retro = True
    else:
        logging.error('Cannot interpret partial_retro: '+str(params.partial_retro))
        exit(1)
    '''
    if not os.path.exists(params.scope_csv):
        logging.error('The scope file cannot be found: '+str(params.scope_csv))
        exit(1)
    '''
    if not os.path.exists(params.rulesfile):
        logging.error('The rules file cannot be found: '+str(params.rulesfile))
        exit(1)
    if not os.path.exists(params.sinkfile):
        logging.error('The sink file cannot be found: '+str(params.sinkfile))
        exit(1)
    ########## handle the call ###########
    with tempfile.TemporaryDirectory() as tmpInputFolder:
        if params.rulesfile_format=='csv':
            logging.debug('Rules file: '+str(params.rulesfile))
            rulesfile = tmpInputFolder+'/rules.csv'
            shutil.copy(params.rulesfile, rulesfile)
            logging.debug('Rules file: '+str(rulesfile))
        elif params.rulesfile_format=='tar':
            with tarfile.open(params.rulesfile) as rf:
                rf.extractall(tmpInputFolder)
            out_file = glob.glob(tmpInputFolder+'/*.csv')
            if len(out_file)>1:
                logging.error('Cannot detect file: '+str(glob.glob(tmpInputFolder+'/*.csv')))
                exit(1)
            elif len(out_file)==0:
                logging.error('The rules tar input is empty')
                exit(1)
            rulesfile = out_file[0]
        else:
            logging.error('Cannot detect the rules_format: '+str(params.rulesfile_format))
            exit(1)
        sourcefile = tmpInputFolder+'/source.csv'
        shutil.copy(params.sourcefile, sourcefile)
        sinkfile = tmpInputFolder+'/sink.csv'
        shutil.copy(params.sinkfile, sinkfile)
        result = run_rp2(sourcefile,
                         sinkfile,
                         rulesfile,
                         params.max_steps,
                         params.topx,
                         params.dmin,
                         params.dmax,
                         params.mwmax_source,
                         params.mwmax_cof,
                         params.timeout,
                         partial_retro)
        ###########################
        if result[1]==b'timeouterror' or result[1]==b'timeoutwarning':
            if not partial_retro:
                logging.error('Timeout of RetroPath2.0 -- Try increasing the timeout limit of the tool')
            else:
                if result[0]==b'':
                    logging.error('Timeout caused RetroPath2.0 to not find any solutions')
                    exit(1)
                else:
                    logging.warning('Timeout of RetroPath2.0 -- Try increasing the timeout limit of the tool')
                    logging.warning('Returning partial results')
                    with open(params.scope_csv, 'wb') as scope_csv:
                        scope_csv.write(result[0])
                    exit(0)
        elif result[1]==b'memwarning' or result[1]==b'memerror':
            if not partial_retro:
                logging.error('RetroPath2.0 has exceeded its memory limit')
                exit(0)
            else:
                if result[0]==b'':
                    logging.error('Memory limit reached by RetroPath2.0 caused it to not find any solutions')
                    exit[0]
                else:
                    logging.warning('RetroPath2.0 has exceeded its memory limit')
                    logging.warning('Returning partial results')
                    with open(params.scope_csv, 'wb') as scope_csv:
                        scope_csv.write(result[0])
                    exit(0)
        elif result[1]==b'sourceinsinkerror':
            logging.error('Source exists in the sink')
            exit(1)
        elif result[1]==b'sourceinsinknotfounderror':
            logging.error('Cannot find the sink-in-source file')
            exit(1)
        elif result[1]==b'ramerror' or result[1]==b'ramwarning':
            logging.error('Memory allocation error')
            exit(1)
        elif result[1]==b'oserror' or result[1]==b'oswarning':
            logging.error('RetroPath2.0 has generated an OS error')
            exit(1)
        elif result[1]==b'noresultwarning':
            if partial_retro:
                if result[0]==b'':
                    logging.error('No results warning caused it to return no results')
                    exit(1)
                else:
                    logging.warning('RetroPath2.0 did not complete successfully')
                    logging.warning('Returning partial results')
                    with open(params.scope_csv, 'wb') as scope_csv:
                        scope_csv.write(result[0])
            else:
                logging.error('RetroPath2.0 did not complete successfully')
                exit(1)
        elif result[1]==b'noresulterror':
            logging.error('Empty results')
            exit(1)
        elif result[1]==b'noerror':
            logging.info('Successfull execution')
            with open(params.scope_csv, 'wb') as scope_csv:
                scope_csv.write(result[0])
        else:
            logging.error('Could not recognise the status message returned: '+str(results[1]))
            exit(1)
