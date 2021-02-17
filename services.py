#!/usr/bin/env python3

"""
Created on December 7 2020

@author: Melchior du Lac
@description:

"""
import json
import os
import tempfile
import logging
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort, Response, make_response
from flask_restful import Resource, Api

from metaxime import rpCache
from metaxime import rpGraph
from metaxime import rpReader
from metaxime import rpFBA
from metaxime import rpSBML
from metaxime import rpEquilibrator
from metaxime import rpSelenzyme
from metaxime import rpGlobalScore

#from logsetup import LOGGING_CONFIG
#logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)

app = Flask(__name__)

GLOBAL_RPCACHE = rpCache()
GLOBAL_RPCACHE.populateCache()
SELENZYME_CACHE_PATH = '/home/metaxime/input_cache/rpselenzyme_data.tar.xz'
SELENZYME_OBJ = rpSelenzyme()
SELENZYME_PC, SELENZYME_UNIPROT_AA_LENGTH, SELENZYNE_DATA_DIR = SELENZYME_OBJ.loadCache(SELENZYME_CACHE_PATH)

#######################################################
##################### HELPER ##########################
#######################################################


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


#######################################################
############## api ###################################
#######################################################


@app.route("/api", methods=["GET", "POST"])
def RestApp():
    """The Flask methods that we support, post and get
    """
    return jsonify(stamp(None))


@app.route("/api/rpReader", methods=["POST"])
def rpReaderService():
    with tempfile.TemporaryDirectory() as tmpdir:
        rp2_file = os.path.join(tmpdir, 'rp2.csv')
        try:
            with open(rp2_file, 'wb') as fo:
                fo.write(request.files['rp2_file'].read())
            rp2paths_compounds_file = os.path.join(tmpdir, 'rp2paths_compounds.csv')
            with open(rp2paths_compounds_file, 'wb') as fo:
                fo.write(request.files['rp2paths_compounds_file'].read())
            rp2paths_pathways_file = os.path.join(tmpdir, 'rp2paths_pathways.csv')
            with open(rp2paths_pathways_file, 'wb') as fo:
                fo.write(request.files['rp2paths_pathways_file'].read())
        except KeyError as e:
            return Response('A required file is missing: '+str(e), status=400)
        rpcollection_file = os.path.join(tmpdir, 'rpcollection.tar.xz')
        status = rpReader.rp2ToCollection(rp2_file,
                                          rp2paths_compounds_file,
                                          rp2paths_pathways_file,
                                          rpcollection_file,
                                          rpcache=GLOBAL_RPCACHE)
        if not status:
            return Response('rpReader has encountered a problem.... please investigate futher', status=400)
        #rpcollection_file.seek(0)
        return send_file(rpcollection_file,
                         as_attachment=True,
                         attachment_filename='rpcollection.tar.xz',
                         mimetype='application/x-tar')


@app.route("/api/rpEquilibrator", methods=["POST"])
def rpEquilibratorService():
    with tempfile.TemporaryDirectory() as tmpdir:
        ############### file in ##############
        rpcollection_file = os.path.join(tmpdir, 'rpcollection.tar.xz')
        try:
            with open(rpcollection_file, 'wb') as fo:
                fo.write(request.files['rpcollection_file'].read())
        except KeyError as e:
            return Response('A required file is missing: '+str(e), status=400)
        ############### parameters #############
        try:
            params = json.load(request.files['data'])
        except ValueError as e:
            return Response('One or more parameters are malformed: '+str(e), status=400)
        except KeyError as e:
            return Response('One or more of the parameters are missing: '+str(e), status=400)
        try:
            ph = float(params['ph'])
        except KeyError:
            app.logger.info('No ph passed. Setting to default to 7.5')
            ph = 7.5
        except ValueError:
            app.logger.warning('ph is not recognised. Setting to default 7.5')
            ph = 7.5
        try:
            temp_k = float(params['temp_k'])
        except KeyError:
            app.logger.info('No temp_k passed. Setting to default to 298.15')
            temp_k = 298.15
        except ValueError:
            app.logger.warning('temp_k is not recognised. Setting to default 298.15')
            temp_k = 298.15
        try:
            ionic_strength = float(params['ionic_strength'])
        except KeyError:
            app.logger.info('No ionic_strength passed. Setting to default to 200.0')
            ionic_strength = 200.0
        except ValueError:
            app.logger.warning('ionic_strength is not recognised. Setting to default 200.0')
            ionic_strength = 200.0
        status = rpEquilibrator.runCollection(rpcollection_file,
                                              rpcollection_file,
                                              ph=ph,
                                              ionic_strength=ionic_strength,
                                              temp_k=temp_k,
                                              rpcache=GLOBAL_RPCACHE)
        if not status:
            return Response('rpEquilibrator has encountered a problem.... please investigate futher', status=400)
        #rpcollection_file.seek(0)
        return send_file(rpcollection_file,
                         as_attachment=True,
                         attachment_filename='rpcollection.tar.xz',
                         mimetype='application/x-tar')


@app.route("/api/rpFBA", methods=["POST"])
def rpFBAService():
    with tempfile.TemporaryDirectory() as tmpdir:
        ########## file in ###############
        rpcollection_file = os.path.join(tmpdir, 'rpcollection.tar.xz')
        gem_file = os.path.join(tmpdir, 'gem_file.sbml')
        try:
            with open(rpcollection_file, 'wb') as fo:
                fo.write(request.files['rpcollection_file'].read())
            with open(gem_file, 'wb') as fo:
                fo.write(request.files['gem_file'].read())
        except KeyError as e:
            return Response('A required file is missing: '+str(e), status=400)
        ########## parameters ############
        try:
            params = json.load(request.files['data'])
        except ValueError as e:
            return Response('One or more parameters are malformed: '+str(e), status=400)
        except KeyError as e:
            return Response('One or more of the parameters are missing: '+str(e), status=400)
        try:
            num_workers = int(params['num_workers'])
        except KeyError:
            app.logger.info('No num_workers passed. Setting to default to 1')
            num_workers = 1
        except ValueError:
            app.logger.warning('num_workers is not recognised. Setting to default 1')
            num_workers = 1
        try:
            keep_merged = bool(params['keep_merged'])
        except KeyError:
            app.logger.info('No keep_merged passed. Setting to default to True')
            keep_merged = True
        except ValueError:
            app.logger.warning('keep_merged is not recognised. Setting to default True')
            keep_merged = True
        try:
            del_sp_pro = bool(params['del_sp_pro'])
        except KeyError:
            app.logger.info('No del_sp_pro passed. Setting to default to False')
            del_sp_pro = False
        except ValueError:
            app.logger.warning('del_sp_pro is not recognised. Setting to default False')
            del_sp_pro = False
        try:
            del_sp_react = bool(params['del_sp_react'])
        except KeyError:
            app.logger.info('No del_sp_react passed. Setting to default to False')
            del_sp_react = False
        except ValueError:
            app.logger.warning('del_sp_react is not recognised. Setting to default False')
            del_sp_react = False
        status = rpFBA.runCollection(rpcollection_file,
                                     gem_file,
                                     rpcollection_file,
                                     num_workers=num_workers,
                                     keep_merged=keep_merged,
                                     del_sp_pro=del_sp_pro,
                                     del_sp_react=del_sp_react,
                                     rpcache=GLOBAL_RPCACHE)
        if not status:
            return Response('rpFBA has encountered a problem.... please investigate futher', status=400)
        #rpcollection_file.seek(0)
        return send_file(rpcollection_file,
                         as_attachment=True,
                         attachment_filename='rpcollection.tar.xz',
                         mimetype='application/x-tar')


@app.route("/api/rpSelenzyme", methods=["POST"])
def rpSelenzymeService():
    with tempfile.TemporaryDirectory() as tmpdir:
        ########### file in ####################
        rpcollection_file = os.path.join(tmpdir, 'rpcollection.tar.xz')
        try:
            with open(rpcollection_file, 'wb') as fo:
                fo.write(request.files['rpcollection_file'].read())
        except KeyError as e:
            return Response('A required file is missing: '+str(e), status=400)
        ########## Parameters ###################
        try:
            params = json.load(request.files['data'])
        except ValueError as e:
            return Response('One or more parameters are malformed: '+str(e), status=400)
        except KeyError as e:
            return Response('One or more of the parameters are missing: '+str(e), status=400)
        try:
            taxo_id = int(params['taxo_id'])
        except KeyError:
            app.logger.info('No taxo_id passed. Setting to default to E.Coli (83333)')
            taxo_id = 83333
        except ValueError:
            app.logger.warning('taxo_id is not recognised. Setting to default E.Coli (83333)')
            taxo_id = 83333
        status = rpSelenzyme.runCollection(rpcollection_file,
                                           taxo_id,
                                           rpcollection_file,
                                           is_cleanup=False,
                                           uniprot_aa_length=SELENZYME_UNIPROT_AA_LENGTH,
                                           data_dir=SELENZYNE_DATA_DIR,
                                           pc=SELENZYME_PC,
                                           rpcache=GLOBAL_RPCACHE)
        #rpcollection_file.seek(0)
        if not status:
            return Response('rpSelenzyme has encountered a problem.... please investigate futher', status=400)
        return send_file(rpcollection_file,
                         as_attachment=True,
                         attachment_filename='rpcollection.tar.xz',
                         mimetype='application/x-tar')


@app.route("/api/rpGlobalScore", methods=["POST"])
def rpGlobalScoreService():
    with tempfile.TemporaryDirectory() as tmpdir:
        rpcollection_file = os.path.join(tmpdir, 'rpcollection.tar.xz')
        try:
            with open(rpcollection_file, 'wb') as fo:
                fo.write(request.files['rpcollection_file'].read())
        except KeyError as e:
            return Response('A required file is missing: '+str(e), status=400)
        status = rpGlobalScore.runCollection(rpcollection_file,
                                             rpcollection_file,
                                             rpcache=GLOBAL_RPCACHE)
        if not status:
            return Response('rpGlobalScore has encountered a problem.... please investigate futher', status=400)
        #rpcollection_file.seek(0)
        return send_file(rpcollection_file,
                         as_attachment=True,
                         attachment_filename='rpcollection.tar.xz',
                         mimetype='application/x-tar')


@app.route("/api/rpPipeline", methods=["POST"])
def rpPipelineService():
    with tempfile.TemporaryDirectory() as tmpdir:
        ########### file in ################
        rp2_file = os.path.join(tmpdir, 'rp2_file.csv')
        rp2paths_compounds_file = os.path.join(tmpdir, 'rp2paths_compounds_file.csv')
        rp2paths_pathways_file = os.path.join(tmpdir, 'rp2paths_pathways_file')
        gem_file = os.path.join(tmpdir, 'gem_file.sbml')
        try:
            with open(rp2_file, 'wb') as fo:
                fo.write(request.files['rp2_file'].read())
            with open(rp2paths_compounds_file, 'wb') as fo:
                fo.write(request.files['rp2paths_compounds_file'].read())
            with open(rp2paths_pathways_file, 'wb') as fo:
                fo.write(request.files['rp2paths_pathways_file'].read())
            with open(gem_file, 'wb') as fo:
                fo.write(request.files['gem_file'].read())
        except KeyError as e:
            return Response('A required file is missing: '+str(e), status=400)
        ############ parameters ##########
        try:
            params = json.load(request.files['data'])
        except ValueError as e:
            return Response('One or more parameters are malformed: '+str(e), status=400)
        except KeyError as e:
            return Response('One or more of the parameters are missing: '+str(e), status=400)
        try:
            ph = float(params['ph'])
        except KeyError:
            app.logger.info('No ph passed. Setting to default to 7.5')
            ph = 7.5
        except ValueError:
            app.logger.warning('ph is not recognised. Setting to default 7.5')
            ph = 7.5
        try:
            temp_k = float(params['temp_k'])
        except KeyError:
            app.logger.info('No temp_k passed. Setting to default to 298.15')
            temp_k = 298.15
        except ValueError:
            app.logger.warning('temp_k is not recognised. Setting to default 298.15')
            temp_k = 298.15
        try:
            ionic_strength = float(params['ionic_strength'])
        except KeyError:
            app.logger.info('No ionic_strength passed. Setting to default to 200.0')
            ionic_strength = 200.0
        except ValueError:
            app.logger.warning('ionic_strength is not recognised. Setting to default 200.0')
            ionic_strength = 200.0
        try:
            num_workers = int(params['num_workers'])
        except KeyError:
            app.logger.info('No num_workers passed. Setting to default to 1')
            num_workers = 1
        except ValueError:
            app.logger.warning('num_workers is not recognised. Setting to default 1')
            num_workers = 1
        try:
            keep_merged = bool(params['keep_merged'])
        except KeyError:
            app.logger.info('No keep_merged passed. Setting to default to True')
            keep_merged = True
        except ValueError:
            app.logger.warning('keep_merged is not recognised. Setting to default True')
            keep_merged = True
        try:
            del_sp_pro = bool(params['del_sp_pro'])
        except KeyError:
            app.logger.info('No del_sp_pro passed. Setting to default to False')
            del_sp_pro = False
        except ValueError:
            app.logger.warning('del_sp_pro is not recognised. Setting to default False')
            del_sp_pro = False
        try:
            del_sp_react = bool(params['del_sp_react'])
        except KeyError:
            app.logger.info('No del_sp_react passed. Setting to default to False')
            del_sp_react = False
        except ValueError:
            app.logger.warning('del_sp_react is not recognised. Setting to default False')
            del_sp_react = False
        try:
            taxo_id = int(params['taxo_id'])
        except KeyError:
            app.logger.info('No taxo_id passed. Setting to default to None to try and recover from GEM SBML')
            taxo_id = None
        except ValueError:
            app.logger.warning('taxo_id is not recognised. Setting to default None to try to recover from GEM SBML')
            taxo_id = None
        rpcollection_file = os.path.join(tmpdir, 'rpcollection.tar.xz')
        rpre_status = rpReader.rp2ToCollection(rp2_file,
                                               rp2paths_compounds_file,
                                               rp2paths_pathways_file,
                                               rpcollection_file,
                                               rpcache=GLOBAL_RPCACHE)
        if not rpre_status:
            return Response('rpReader has encountered a problem.... please investigate futher', status=400)
        rpeq_status = rpEquilibrator.runCollection(rpcollection_file,
                                                   rpcollection_file,
                                                   ph=ph,
                                                   ionic_strength=ionic_strength,
                                                   temp_k=temp_k,
                                                   rpcache=GLOBAL_RPCACHE)
        if not rpeq_status:
            return Response('rpEquilibrator has encountered a problem.... please investigate futher', status=400)
        rpfba_status = rpFBA.runCollection(rpcollection_file,
                                           gem_file,
                                           rpcollection_file,
                                           num_workers=num_workers,
                                           keep_merged=keep_merged,
                                           del_sp_pro=del_sp_pro,
                                           del_sp_react=del_sp_react,
                                           rpcache=GLOBAL_RPCACHE)
        if not rpfba_status:
            return Response('rpFBA has encountered a problem.... please investigate futher', status=400)
        if taxo_id==None:
            #if you cannot find the annotation then try to recover it from the GEM file
            rpsbml_gem = rpSBML(model_name='tmp', path=gem_file)
            taxo_id = rpsbml_gem.readTaxonomy() #TODO fix, seems to be a list
            logging.info('The taxonomy_id is '+str(taxo_id))
        rpsel_status = rpSelenzyme.runCollection(rpcollection_file,
                                                 taxo_id,
                                                 rpcollection_file,
                                                 is_cleanup=False,
                                                 uniprot_aa_length=SELENZYME_UNIPROT_AA_LENGTH,
                                                 data_dir=SELENZYNE_DATA_DIR,
                                                 pc=SELENZYME_PC,
                                                 rpcache=GLOBAL_RPCACHE)
        if not rpsel_status:
            return Response('rpSelenzyme has encountered a problem.... please investigate futher', status=400)
        rpglo_status = rpGlobalScore.runCollection(rpcollection_file,
                                                   rpcollection_file,
                                                   rpcache=GLOBAL_RPCACHE)
        if not rpglo_status:
            return Response('rpGlobalScore has encountered a problem.... please investigate futher', status=400)
        #rpcollection_file.seek(0)
        return send_file(rpcollection_file,
                         as_attachment=True,
                         attachment_filename='rpcollection.tar.xz',
                         mimetype='application/x-tar')

if __name__== "__main__":
    #handler = RotatingFileHandler('metaxime.log', maxBytes=10000, backupCount=1)
    #handler.setLevel(logging.DEBUG)
    #logger.addHandler(handler)
    app.run(host="0.0.0.0", port=7777, debug=True, threaded=False)
