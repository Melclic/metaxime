#!/usr/bin/env python3

"""
Created on March 7 2019

@author: Melchior du Lac
@description: Standalone version of RP2paths. Returns bytes to be able to use the same file in REST application

"""
import subprocess
import resource
import tempfile
import glob
import io
import logging
import argparse


logging.basicConfig(
    #level=logging.DEBUG,
    #level=logging.WARNING,
    level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

MAX_VIRTUAL_MEMORY = 20000 * 1024 * 1024 # 20GB -- define what is the best
#MAX_VIRTUAL_MEMORY = 20 * 1024 * 1024 # 20GB -- define what is the best

##
#
#
def limit_virtual_memory():
    """The function to set the memory limits
    """
    resource.setrlimit(resource.RLIMIT_AS, (MAX_VIRTUAL_MEMORY, resource.RLIM_INFINITY))

def run_rp2paths(rp2_pathways, timeout, logger=None):
    """Call the KNIME RetroPath2.0 workflow

    :param rp2_pathways: The path to the RetroPath2.0 scope results
    :param timeout: The timeout of the function in minutes
    :param logger: Logger object (Default: None)

    :param source_bytes: str
    :param sink_bytes: int
    :param logger: logging

    :rtype: tuple
    :return: tuple of bytes with the out_paths results, compounds results, the status message, the command used
    """
    ### not sure why throws an error:
    if logger==None:
        logging.basicConfig(level=logging.DEBUG)
        logger = logging.getLogger(__name__)
    out_paths = b''
    out_compounds = b''
    with tempfile.TemporaryDirectory() as tmp_output_folder:
        rp2paths_command = 'python /home/RP2paths.py all '+str(rp2_pathways)+' --outdir '+str(tmp_output_folder)+' --timeout '+str(int(timeout*60.0))
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
                return b'', b'', b'timeout', str.encode('Command: '+str(rp2paths_command)+'\n Error: '+str(error)+'\n tmp_output_folder: '+str(glob.glob(os.path.join(tmp_output_folder, '*'))))
            if 'failed to map segment from shared object' in error:
                logger.error('RP2paths does not have sufficient memory to continue')
                return b'', b'', b'memoryerror', str.encode('Command: '+str(rp2paths_command)+'\n Error: '+str(error)+'\n tmp_output_folder: '+str(glob.glob(os.path.join(tmp_output_folder, '*'))))
            ### convert the result to binary and return ###
            try:
                with open(os.path.join(tmp_output_folder, 'out_paths.csv'), 'rb') as op:
                    out_paths = op.read()
                with open(os.path.join(tmp_output_folder, 'compounds.txt'), 'rb') as c:
                    out_compounds = c.read()
                return out_paths, out_compounds, b'noerror', b''
            except FileNotFoundError as e:
                logger.error('Cannot find the output files out_paths.csv or compounds.txt')
                return b'', b'', b'filenotfounderror', str.encode('Command: '+str(rp2paths_command)+'\n Error: '+str(e)+'\n tmp_output_folder: '+str(glob.glob(os.path.join(tmp_output_folder, '*'))))
        except OSError as e:
            logger.error('Subprocess detected an error when calling the rp2paths command')
            return b'', b'', b'oserror', str.encode('Command: '+str(rp2paths_command)+'\n Error: '+str(e)+'\n tmp_output_folder: '+str(glob.glob(os.path.join(tmp_output_folder, '*'))))
        except ValueError as e:
            logger.error('Cannot set the RAM usage limit')
            return b'', b'', b'ramerror', str.encode('Command: '+str(rp2paths_command)+'\n Error: '+str(e)+'\n tmp_output_folder: '+str(glob.glob(os.path.join(tmp_output_folder, '*'))))

# Wrapper for the RP2paths script that takes the same input (results.csv) as the original script but returns
# the out_paths.csv so as to be compliant with Galaxy
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wrapper for the python RP2paths script')
    parser.add_argument('-rp_pathways', type=str)
    parser.add_argument('-rp2paths_pathways', type=str)
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-timeout', type=int, default=30)
    params = parser.parse_args()
    if params.timeout<=0:
        logging.error('Timeout cannot be less or equal to 0 :'+str(params.timeout))
        exit(1)
    result = run_rp2paths(params.rp_pathways, params.timeout)
    if result[2]==b'timeout':
        logging.error('RP2paths has reached its timeout limit')
        exit(1)
    if result[2]==b'filenotfounderror':
        logging.error('Cannot detect the out_path.csv and/or compounds.txt files')
        exit(1)
    elif result[2]==b'memoryerror':
        logging.error('Memory allocation error')
        exit(1)
    elif result[2]==b'oserror':
        logging.error('rp2paths has generated an OS error')
        exit(1)
    elif result[2]==b'ramerror':
        logging.error('Could not setup a RAM limit')
        exit(1)
    elif result[0]==b'':
        logging.error('Empty rp2paths_pathways')
        exit(1)
    elif result[1]==b'':
        logging.error('Empty rp2paths_compounds')
        exit(1)
    with open(params.rp2paths_pathways, 'wb') as out_paths:
        out_paths.write(result[0])
    with open(params.rp2paths_compounds, 'wb') as out_compounds:
        out_compounds.write(result[1])
