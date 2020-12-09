#!/usr/bin/env python3

"""
Created on December 7 2020

@author: Melchior du Lac
@description:

"""
import json
import os
import csv
import tempfile
import time
import sys

import logging
'''
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort, Response, make_response
from flask_restful import Resource, Api
from logging.handlers import RotatingFileHandler

from logging.config import dictConfig

dictConfig({
    'version': 0.1,
    'formatters': {'default': {
        'format': '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
    }},
    'handlers': {'wsgi': {
        'class': 'logging.StreamHandler',
        'stream': 'ext://flask.logging.wsgi_errors_stream',
        'formatter': 'default'
    }},
    'root': {
        'level': 'DEBUG',
        'handlers': ['wsgi']
    }
})


from rq import Connection, Queue
from redis import Redis


#######################################################
############## REST ###################################
#######################################################

app = Flask(__name__)
api = Api(app)

#app.logger.setLevel(logging.WARNING)

def stamp(data, status=1):
    """Default message to return

    :param data: The data to be passes
    :param status: The int value of the status

    :type data: dict
    :type status: int

    :rtype: dict
    :return: The dict of the stamp
    """
    appinfo = {'app': 'Metaxime', 'version': '0.1',
               'author': 'Melchior du Lac',
               'time': datetime.now().isoformat(),
               'status': status}
    out = appinfo.copy()
    out['data'] = data
    return out


class RestApp(Resource):
    """The Flask methods that we support, post and get
    """
    def post(self):
        return jsonify(stamp(None))
    def get(self):
        return jsonify(stamp(None))


# NOTE: Avoid returning numpy or pandas object in order to keep the client lighter.
class RestQuery(Resource):
    """Class containing the REST requests for RP2
    """
    def post(self):
        """Make the REST request using the POST method

        :rtype: Response
        :return: Flask Response object
        """
        source_file_bytes = request.files['source_file'].read()
        sink_file_bytes = request.files['sink_file'].read()
        rules_file_bytes = request.files['rules_file'].read()
        params = json.load(request.files['data'])
        ##### REDIS ##############
        conn = Redis()
        q = Queue('default', connection=conn, default_time_out='24h')
        #pass the cache parameters to the rpCofactors object
        if params['partial_retro']=='True' or params['partial_retro']=='T' or params['partial_retro']=='true' or params['partial_retro']==True:
            partial_retro = True
        elif params['partial_retro']=='True' or params['partial_retro']=='F' or params['partial_retro']=='false' or params['partial_retro']==False:
            partial_retro = False
        else:
            app.logger.warning('Cannot interpret partial_retro: '+str(params['partial_retro']))
            app.logger.warning('Setting to False')
            partial_retro = False
        app.logger.debug('max_steps: '+str(params['max_steps']))
        app.logger.debug('topx: '+str(params['topx']))
        app.logger.debug('dmin: '+str(params['dmin']))
        app.logger.debug('dmax: '+str(params['dmax']))
        app.logger.debug('mwmax_source: '+str(params['mwmax_source']))
        app.logger.debug('mwmax_cof: '+str(params['mwmax_cof']))
        app.logger.debug('time_out: '+str(params['time_out']))
        app.logger.debug('ram_limit: '+str(params['ram_limit']))
        app.logger.debug('partial_retro: '+str(params['partial_retro']))
        async_results = q.enqueue(rpTool.run_rp2,
                                  source_file_bytes,
                                  sink_file_bytes,
                                  rules_file_bytes,
                                  int(params['max_steps']),
                                  int(params['topx']),
                                  int(params['dmin']),
                                  int(params['dmax']),
                                  int(params['mwmax_source']),
                                  int(params['mwmax_cof']),
                                  int(params['time_out']),
                                  int(params['ram_limit']),
                                  partial_retro)
        result = None
        while result is None:
            result = async_results.return_value
            if async_results.get_status()=='failed':
                return Response('Redis job failed \n '+str(result), status=500)
            time.sleep(2.0)
        ###########################
        if result[1]=='timeouterror' or result[1]=='timeoutwarning':
            #for debugging
            app.logger.warning(result[2])
            if not partial_retro:
                app.logger.error('Timeout of RetroPath2.0 -- Try increasing the time_out limit of the tool')
                return Response('Timeout of RetroPath2.0--Try increasing the time_out limit of the tool', status=408)
            else:
                if result[0]=='':
                    return Response('Timeout caused RetroPath2.0 to not find any solutions', status=404)
                else:
                    app.logger.warning('Timeout of RetroPath2.0 -- Try increasing the time_out limit of the tool')
                    app.logger.warning('Returning partial results')
                    status_message = 'WARNING: Timeout of RetroPath2.0--Try increasing the time_out limit of the tool--Returning partial results'
                    scope_csv = io.BytesIO()
                    scope_csv.write(result[0])
                    ###### IMPORTANT ######
                    scope_csv.seek(0)
                    #######################
                    response = make_response(send_file(scope_csv, as_attachment=True, attachment_filename='rp2_pathways.csv', mimetype='text/csv'))
                    response.headers['status_message'] = status_message
                    return response
        elif result[1]=='memwarning' or result[1]=='memerror':
            #for debugging
            app.logger.warning(result[2])
            if not partial_retro:
                app.logger.error('RetroPath2.0 has exceeded its memory limit')
                return Response('RetroPath2.0 has exceeded its memory limit', status=403)
            else:
                if result[0]=='':
                    return Response('Memory limit reached by RetroPath2.0 caused it to not find any solutions', status=404)
                else:
                    app.logger.warning('RetroPath2.0 has exceeded its memory limit')
                    app.logger.warning('Returning partial results')
                    status_message = 'WARNING: RetroPath2.0 has exceeded its memory limit--Returning partial results'
                    scope_csv = io.BytesIO()
                    scope_csv.write(result[0])
                    ###### IMPORTANT ######
                    scope_csv.seek(0)
                    #######################
                    response = make_response(send_file(scope_csv, as_attachment=True, attachment_filename='rp2_pathways.csv', mimetype='text/csv'))
                    response.headers['status_message'] = status_message
                    return response
        elif result[1]=='sourceinsinkerror':
            app.logger.error('Source exists in the sink')
            return Response('Source exists in the sink', status=403)
        elif result[1]=='sourceinsinknotfounderror':
            app.logger.error('Cannot find the sink-in-source file')
            return Response('Cannot find the sink-in-source file', status=500)
        elif result[1]=='ramerror' or result[1]=='ramwarning':
            app.logger.error('Memory allocation error')
            return Response('Memory allocation error', status=500)
        elif result[1]=='oserror' or result[1]=='oswarning':
            app.logger.error('RetroPath2.0 has generated an OS error')
            return Response('RetroPath2.0 returned an OS error', status=500)
        elif result[1]=='noresultwarning':
            if partial_retro:
                if result[0]=='':
                    return Response('No results warning caused it to return no results', status=404)
                else:
                    app.logger.warning('RetroPath2.0 did not complete successfully')
                    app.logger.warning('Returning partial results')
                    status_message = 'WARNING: RetroPath2.0 did not complete successfully--Returning partial results'
                    scope_csv = io.BytesIO()
                    scope_csv.write(result[0])
                    ###### IMPORTANT ######
                    scope_csv.seek(0)
                    #######################
                    response = make_response(send_file(scope_csv, as_attachment=True, attachment_filename='rp2_pathways.csv', mimetype='text/csv'))
                    response.headers['status_message'] = status_message
                    return response
            else:
                return Response('RetroPath2.0 could not complete successfully', status=404)
        elif result[1]=='noresulterror':
            app.logger.error('Empty results')
            return Response('RetroPath2.0 cannot not find any solutions--Try reducing the complexity of the problem', status=404)
        elif result[1]=='noerror':
            status_message = 'Successfull execution'
            scope_csv = io.BytesIO()
            scope_csv.write(result[0])
            ###### IMPORTANT ######
            scope_csv.seek(0)
            #######################
            response = make_response(send_file(scope_csv, as_attachment=True, attachment_filename='rp2_pathways.csv', mimetype='text/csv'))
            response.headers['status_message'] = status_message
            return response
        else:
            app.logger.error('Could not recognise the status message returned: '+str(results[1]))
            return Response('Could not recognise the status message returned: '+str(results[1]), status=500)


api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')

'''
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
             dmin=0,
             dmax=1000,
             mwmax_source=1000,
             mwmax_cof=1000,
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
            return False
        sink_file = os.path.join(tmp_output_folder, 'sink.csv')
        gensink_status = rpExtractSink.genSink(gem_file, sink_file, remove_dead_end=True)
        if not gensink_status:
            return False
        ########## reaction rules ##########
        rules_file = os.path.join(tmp_output_folder, 'reaction_rules.csv')
        rr_status = callRR.passRules(rules_file, rules_type, rules_diameters, 'csv')
        if not rr_status:
            return False
        ########## Retropath2 ##############
        rp2_file = os.path.join(tmp_output_folder, 'rp2_pathways.csv')
        rp2_status = callRP2.run(rp2_file,
                                 source_file,
                                 sink_file,
                                 rules_file,
                                 max_steps,
                                 topx,
                                 dmin,
                                 dmax,
                                 mwmax_source,
                                 mwmax_cof,
                                 timeout,
                                 partial_retro)
        if not rp2_status:
            return False
        ######### rp2paths ################
        rp2paths_pathways_file = os.path.join(tmp_output_folder, 'rp2paths_pathways.csv')
        rp2paths_compounds_file = os.path.join(tmp_output_folder, 'rp2paths_compounds.tsv')
        rp2paths_status = callRP2paths.run(rp2_file, rp2paths_pathways_file, rp2paths_compounds_file)
        if not rp2paths_status:
            return False
        ####### Analysis #################
        rpReader.rp2ToCollection(rp2_file,
                                 rp2paths_compounds_file,
                                 rp2paths_pathways_file,
                                 rpcollection_file,
                                 rpcache=global_rpcache)
        rpFBA.runCollection(rpcollection_file,
                            gem_file,
                            rpcollection_file,
                            num_workers=1,
                            keep_merged=True,
                            del_sp_pro=False,
                            del_sp_react=False,
                            rpcache=global_rpcache)
        rpEquilibrator.runCollection(rpcollection_file,
                                     rpcollection_file,
                                     cc=global_cc,
                                     ph=ph,
                                     ionic_strength=ionic_strength,
                                     temp_k=temp_k,
                                     rpcache=global_rpcache)
        rpSelenzyme.runCollection(rpcollection_file,
                                  taxo_id,
                                  rpcollection_file,    
                                  cache_path=os.path.join('/home', 'metaxime', 'input_cache', 'rpselenzyme_data.tar.xz'),
                                  rpcache=global_rpcache)
        rpGlobalScore.runCollection(rpcollection_file,
                                    rpcollection_file,
                                    rpcache=global_rpcache)

"""
if __name__== "__main__":
    handler = RotatingFileHandler('metaxime.log', maxBytes=10000, backupCount=1)
    handler.setLevel(logging.ERROR)
    app.logger.addHandler(handler)
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=True)
"""
