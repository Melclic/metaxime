#!/usr/bin/env python3
"""
Created on 9 September 2020

@author: Melchior du Lac
@description: Generate SBtab files for Equilibrator

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
         pathway_id='rp_pathway',
         fba_id=None,
         thermo_id='dfG_prime_o',
         ph=7.0,
         ionic_strength=200.0,
         pMg=10.0,
         temp_k=298.15):
    """Given a tar input file, generate the SBtab for each rpSBML file by calling the docker

    :param inputfile: The path to the input file (Either TAR of rpSBML or single rpSBML)
    :param input_format: Either TAR and rpSBML (Valid options: [tar, sbml])
    :param output: The path to the output file
    :param pathway_id: The id of the heterologous pathway of interest (Default: rp_pathway)
    :param fba_id: The id of the FBA objective. If None all fluxes are set to 1.0 (Default: None)
    :param thermo_id: The thermodynamics used (Default: dfG_prime_o)
    :param ph: The pH of the host organism (Default: 7.0)
    :param ionic_strength: Ionic strenght of the host organism (Default: 200.0)
    :param pMg: The pMg of the host organism (Default: 10.0)
    :param temp_k: The temperature of the host organism in Kelvin (Default: 298.15)

    :type inputfile: str
    :type input_format: str
    :type output: str
    :type pathway_id: str
    :type fba_id: str
    :type thermo_id: str
    :type ph: float
    :type ionic_strength: float
    :type pMg: float
    :type temp_k: float

    :rtype: None
    :return: None
    """
    docker_client = docker.from_env()
    image_str = 'brsynth/rpthermo-standalone:v2'
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
        if os.path.exists(inputfile):
            shutil.copy(inputfile, tmpOutputFolder+'/input.dat')
            command = ['/home/tool_eqSBtab.py',
                       '-input',
                       '/home/tmp_output/input.dat',
                       '-input_format',
                       str(input_format),
                       '-pathway_id',
                       str(pathway_id),
                       '-thermo_id',
                       str(thermo_id),
                       '-fba_id',
                       str(fba_id),
                       '-ph',
                       str(ph),
                       '-ionic_strength',
                       str(ionic_strength),
                       '-pMg',
                       str(pMg),
                       '-temp_k',
                       str(temp_k),
                       '-output',
                       '/home/tmp_output/output.dat']
            container = docker_client.containers.run(image_str, 
                                                     command, 
                                                     detach=True, 
                                                     stderr=True, 
                                                     volumes={tmpOutputFolder+'/': {'bind': '/home/tmp_output', 'mode': 'rw'}})
            container.wait()
            err = container.logs(stdout=False, stderr=True)
            err_str = err.decode('utf-8')
            print(err_str)
            if 'ERROR' in err_str:
                print(err_str)
            if 'WARNING' in err_str:
                print(err_str)
            if not os.path.exists(tmpOutputFolder+'/output.dat'):
                logging.error('Did not generate an output file')
            else:
                shutil.copy(tmpOutputFolder+'/output.dat', output)
            container.remove()
        else:
            logging.error('Cannot find the input file: '+str(inputfile))




##
#
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser('Python wapper to execute the MDF on rpSBML heterologous pathways')
    parser.add_argument('-input', type=str)
    parser.add_argument('-output', type=str)
    parser.add_argument('-input_format', type=str)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-fba_id', type=str, default='fba_obj_fraction')
    parser.add_argument('-thermo_id', type=str, default='dfG_prime_o')
    parser.add_argument('-ph', type=float, default=7.0)
    parser.add_argument('-ionic_strength', type=float, default=200.0)
    parser.add_argument('-pMg', type=float, default=10.0)
    parser.add_argument('-temp_k', type=float, default=298.15)
    params = parser.parse_args()
    main(params.input, params.input_format, params.output, params.pathway_id, params.fba_id, params.thermo_id, params.ph, params.ionic_strength, params.pMg, params.temp_k)
