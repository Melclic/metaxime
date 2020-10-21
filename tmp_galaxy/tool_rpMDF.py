#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Galaxy script to query rpThermo docker

"""
import argparse
import shutil
import tarfile
import tempfile
import os
import logging
import sys

logging.basicConfig(
    #level=logging.DEBUG,
    level=logging.WARNING,
    #level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

sys.path.insert(0, '/home/')
import rpToolServe

if __name__ == "__main__":
    parser = argparse.ArgumentParser('Calculate the Min-Max Driving Force (MDF) of heterologous pathways')
    parser.add_argument('-input', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-input_format', type=str)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-fba_id', type=str, default='fba_obj_fraction')
    parser.add_argument('-thermo_id', type=str, default='dfG_prime_o')
    parser.add_argument('-ph', type=float, default=7.5)
    parser.add_argument('-ionic_strength', type=float, default=200.0)
    parser.add_argument('-pMg', type=float, default=10.0)
    parser.add_argument('-temp_k', type=float, default=298.15)
    parser.add_argument('-stdev_factor', type=float, default=1.96)
    params = parser.parse_args()
    if params.input_format=='tar':
        rpToolServe.runMDF_hdd(params.input, params.output, params.pathway_id, params.thermo_id, params.fba_id, params.ph, params.ionic_strength, params.pMg, params.temp_k, params.stdev_factor)
    elif params.input_format=='sbml':
        with tempfile.TemporaryDirectory() as tmpOutputFolder:
            inputTar = tmpOutputFolder+'/tmp_input.tar.xz'
            outputTar = tmpOutputFolder+'/tmp_output.tar.xz'
            with tarfile.open(inputTar, mode='w:xz') as tf:
                info = tarfile.TarInfo('single_rpsbml.xml') #need to change the name since galaxy creates .dat files
                info.size = os.path.getsize(params.input)
                tf.addfile(tarinfo=info, fileobj=open(params.input, 'rb'))
            rpToolServe.runMDF_hdd(inputTar, outputTar, params.pathway_id, params.thermo_id, params.fba_id, params.ph, params.ionic_strength, params.pMg, params.temp_k, params.stdev_factor)
            with tarfile.open(outputTar) as outTar:
                outTar.extractall(tmpOutputFolder)
            out_file = glob.glob(tmpOutputFolder+'/*.xml')
            if len(out_file)>1:
                logging.warning('There are more than one output file...')
            shutil.copy(out_file[0], params.output)
    else:
        logging.error('Cannot identify the input/output format: '+str(params.input_format))
