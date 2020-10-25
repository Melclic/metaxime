#!/usr/bin/env python3
"""
Created on September 21 2019

@author: Melchior du Lac
@description: Extract the sink from an SBML into RP2 friendly format

"""
import argparse
import tempfile
import os
import logging
import shutil
import docker


##
#
#
def main(inputfile,
         input_format,
         output,
         weight_rp_steps=0.10002239003499142,
         weight_rule_score=0.13346271414277305,
         weight_fba=0.6348436269211155,
         weight_thermo=0.13167126890112002,
         max_rp_steps=15,
         topX=10,
         thermo_ceil=5000.0,
         thermo_floor=-5000.0,
         fba_ceil=5.0,
         fba_floor=0.0,
         pathway_id='rp_pathway',
         objective_id='obj_fraction',
         thermo_id='dfG_prime_m'):
    docker_client = docker.from_env()
    image_str = 'brsynth/rpglobalscore-standalone:v2'
    try:
        image = docker_client.images.get(image_str)
    except docker.errors.ImageNotFound:
        logging.warning('Could not find the image, trying to pull it')
        try:
            docker_client.images.pull(image_str)
            image = docker_client.images.get(image_str)
        except docker.errors.ImageNotFound:
            logging.error('Cannot pull image: '+str(image_str))
            exit(1)
    with tempfile.TemporaryDirectory() as tmpOutputFolder:
        shutil.copy(inputfile, tmpOutputFolder+'/input.dat')
        command = ['/home/tool_rpGlobalScore.py',
                   '-input',
                   '/home/tmp_output/input.dat',
                   '-input_format',
                   input_format,
                   '-output',
                   '/home/tmp_output/output.dat',
                   '-weight_rule_score',
                   weight_rule_score,
                   '-weight_fba',
                   weight_fba,
                   '-weight_thermo',
                   weight_thermo,
                   '-weight_rp_steps',
                   weight_rp_steps,
                   '-max_rp_steps',
                   max_rp_steps,
                   '-topX',
                   topX,
                   '-thermo_ceil',
                   thermo_ceil,
                   '-thermo_floor',
                   thermo_floor,
                   '-fba_ceil',
                   fba_ceil,
                   '-fba_floor',
                   fba_floor,
                   '-pathway_id',
                   pathway_id,
                   '-objective_id',
                   objective_id,
                   '-thermo_id',
                   thermo_id]
        command = [str(i) for i in command]
        docker_client.containers.run(image_str, 
                command, 
                auto_remove=True, 
                detach=False, 
                volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
        shutil.copy(tmpOutputFolder+'/output.dat', output)


##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Generate the sink from a model SBML by specifying the compartment')
    parser.add_argument('-input', type=str)
    parser.add_argument('-input_format', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-weight_rule_score', type=float, default=0.13346271414277305)
    parser.add_argument('-weight_fba', type=float, default=0.6348436269211155)
    parser.add_argument('-weight_thermo', type=float, default=0.13167126890112002)
    parser.add_argument('-weight_rp_steps', type=float, default=0.10002239003499142)
    parser.add_argument('-max_rp_steps', type=int, default=15) #WARNING: should not have a default
    parser.add_argument('-topX', type=int, default=10)
    parser.add_argument('-thermo_ceil', type=float, default=5000.0)
    parser.add_argument('-thermo_floor', type=float, default=-5000.0)
    parser.add_argument('-fba_ceil', type=float, default=5.0)
    parser.add_argument('-fba_floor', type=float, default=0.0)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-objective_id', type=str, default='obj_fraction')
    parser.add_argument('-thermo_id', type=str, default='dfG_prime_m')
    params = parser.parse_args()
    main(params.input,
         params.input_format,
         params.output,
         params.weight_rule_score,
         params.weight_fba,
         params.weight_thermo,
         params.weight_rp_steps,
         params.max_rp_steps,
         params.topX,
         params.thermo_ceil,
         params.thermo_floor,
         params.fba_ceil,
         params.fba_floor,
         params.pathway_id,
         params.objective_id,
         params.thermo_id)
