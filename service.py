#!/usr/bin/env python3

"""
Created on December 7 2020

@author: Melchior du Lac
@description:

"""
import json
import os
import glob
import tempfile
import time
import io
import collections.abc

import logging
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort, Response, make_response
from flask_restful import Resource, Api
from logging.handlers import RotatingFileHandler


from rq import Connection, Queue
from redis import Redis
from rq.job import Job

import pipeline
from flask_cors import CORS

from logsetup import LOGGING_CONFIG
logging.config.dictConfig(LOGGING_CONFIG)
logger = logging.getLogger(__name__)

'''
from logging.config import dictConfig
dictConfig({
    'version': 1,
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
'''

app = Flask(__name__)
CORS(app)

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


#######################################################
############## REST ###################################
#######################################################


@app.route("/REST", methods=["GET", "POST"])
def RestApp():
    """The Flask methods that we support, post and get
    """
    return jsonify(stamp(None))


@app.route("/REST/SubmitJob", methods=["GET", "POST"])
def submitJob():
    """Make the REST request using the POST method

    :rtype: Response
    :return: Flask Response object
    """
    #params = json.load(request.files['data'])
    params = request.get_json()
    logger.debug(params)
    logger.debug('name: '+str(params['name']))
    logger.debug('SMILES: '+str(params['smiles']))
    logger.debug('GEM: '+str(params['gem']))
    logger.debug('Steps: '+str(params['steps']))
    if not params['smiles'] or not params['gem'] or not params['steps']:
        response = make_response(jsonify({}), 500)
        response.headers["Content-Type"] = "application/json"
        return response
    ##### REDIS ##############
    redis_conn = Redis()
    q = Queue('default', connection=redis_conn, default_time_out='24h')
    job = q.enqueue(pipeline.pipeline, str(params['smiles']), str(params['gem']), int(params['steps']), job_timeout='24h')
    job.meta['progress'] = 0
    job.meta['err_msg'] = ''
    job.save_meta()
    logger.debug('job: '+str(job))
    logger.debug('job.id: '+str(job.id))
    modResJSON(job.id,
               job.get_status(),
               job.meta,
               params['name'],
               params['smiles'],
               params['gem'],
               params['steps'])
    ####### return info ########
    toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta}
    logger.debug(toret)
    response = make_response(jsonify(toret), 201)
    response.headers["Content-Type"] = "application/json"
    return response


@app.route("/REST/RetreiveJob", methods=["GET", "POST"])
def retreiveJob():
    params = json.load(request.files['data'])
    redis_conn = Redis()
    logger.debug('Retreiving results from job: '+str(params['job_id']))
    job = Job.fetch(params['job_id'], connection=redis_conn)
    logger.debug('job status: '+str(job.get_status()))
    if job.get_status()=='failed':
        logger.error('The job '+str(job.id)+' failed: '+str(job.return_value))
        toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta}
        modResJSON(job.id, job.get_status(), job.meta)
        logger.debug(toret)
        response = make_response(jsonify(toret), 500)
        response.headers["Content-Type"] = "application/json"
        return response
    elif job.get_status()=='finished':
        pipeline_result = job.return_value
        if not pipeline_result[0]:
            logger.error('The job '+str(job.id)+' failed: '+str(pipeline_result[1]))
            toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta, 'pipeline_status': pipeline_result[0], 'pipeline_error_msg': pipeline_result[2]}
            modResJSON(job.id, job.get_status(), job.meta)
            logger.debug(toret)
            response = make_response(jsonify(toret), 500)
            response.headers["Content-Type"] = "application/json"
            return response
        else:
            logger.debug('Successfull run')
            modResJSON(job.id, job.get_status(), job.meta)
            binary_rpcol = io.BytesIO()
            binary_rpcol.write(pipeline_result[2])
            ###### IMPORTANT ######
            binary_rpcol.seek(0)
            #######################
            response = make_response(send_file(binary_rpcol, as_attachment=True, attachment_filename=str(job.id)+'.tar.xz', mimetype='application/x-tar'), 200)
            response.headers['status_message'] = 'success'
            return response
    else:
        logger.debug('The job does not seem to be ready: '+str(job.get_status()))
        toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta}
        modResJSON(job.id, job.get_status(), job.meta)
        logger.debug(toret)
        response = make_response(jsonify(toret), 102)
        response.headers["Content-Type"] = "application/json"
        return response


@app.route("/REST/QueryJob", methods=["GET", "POST"])
def queryJob():
    params = json.load(request.files['data'])
    redis_conn = Redis()
    logger.debug('Received a query for job: '+str(params['job_id']))
    job = Job.fetch(params['job_id'], connection=redis_conn)
    logger.debug(job)
    toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta}
    modResJSON(job.id, job.get_status(), job.meta)
    logger.debug(toret)
    response = make_response(jsonify(toret), 200)
    response.headers["Content-Type"] = "application/json"
    return response
        

@app.route("/REST/GetJobList", methods=["GET", "POST"])
def getJobList():
    params = request.get_json()
    mx_res = {}
    with open('/mx-results/job_results.json', 'r') as jr:
        mx_res = json.load(jr)
    response = make_response(jsonify(mx_res), 200)
    response.headers["Content-Type"] = "application/json"
    return response

@app.route("/REST/GetResults", methods=["GET", "POST"])
def getResults():
    params = request.get_json()
    mx_res = []
    for i in glob.glob(os.path.join('/mx-results/', params['job_id'], 'rpsbml_collection', 'model_json', '*')):
        with open(i, 'r') as j:
            mx_res.append(j.read())
    response = make_response(jsonify(mx_res), 200)
    response.headers["Content-Type"] = "application/json"
    return response

@app.route("/REST/GetNetwork", methods=["GET", "POST"])
def getNetwork():
    params = request.get_json()
    mx_res = {}
    with open(os.path.join('/mx-results/', params['job_id'], 'rpsbml_collection', 'networks', params['rpsbml_id']+'.json'), 'r') as jr:
        mx_res = json.load(jr)
    response = make_response(jsonify(mx_res), 200)
    response.headers["Content-Type"] = "application/json"
    return response


if __name__== "__main__":
    #handler = RotatingFileHandler('metaxime.log', maxBytes=10000, backupCount=1)
    #handler.setLevel(logging.DEBUG)
    #logger.addHandler(handler)
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=False)
