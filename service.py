#!/usr/bin/env python3

"""
Created on December 7 2020

@author: Melchior du Lac
@description:

"""
import json
import os
import tempfile
import time
import io

import logging
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort, Response, make_response
from flask_restful import Resource, Api
from logging.handlers import RotatingFileHandler


from rq import Connection, Queue
from redis import Redis
from rq.job import Job

import pipeline

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

#######################################################
############## REST ###################################
#######################################################

app = Flask(__name__)
api = Api(app)

#logger.setLevel(logging.WARNING)

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
class RestSubmitJob(Resource):
    """Class containing the REST requests for the pipeline
    """
    def post(self):
        """Make the REST request using the POST method

        :rtype: Response
        :return: Flask Response object
        """
        params = json.load(request.files['data'])
        logger.debug('name: '+str(params['name']))
        logger.debug('SMILES: '+str(params['smiles']))
        logger.debug('GEM: '+str(params['gem']))
        logger.debug('Steps: '+str(params['steps']))
        ##### REDIS ##############
        redis_conn = Redis()
        q = Queue('default', connection=redis_conn, default_time_out='24h')
        job = q.enqueue(pipeline.pipeline, str(params['smiles']), str(params['gem']), int(params['steps']), job_timeout='24h')
        job.meta['progress'] = 0
        job.save_meta()
        logger.debug('job: '+str(job))
        logger.debug('job.id: '+str(job.id))
        toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta}
        logger.debug(toret)
        response = make_response(jsonify(toret), 201)
        response.headers["Content-Type"] = "application/json"
        return response


class RestRetreiveJob(Resource):
    def post(self):
        params = json.load(request.files['data'])
        redis_conn = Redis()
        logger.debug('Retreiving results from job: '+str(params['job_id']))
        job = Job.fetch(params['job_id'], connection=redis_conn)
        logger.debug('job status: '+str(job.get_status()))
        if job.get_status()=='failed':
            logger.error('The job '+str(job.id)+' failed: '+str(job.return_value))
            toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta}
            logger.debug(toret)
            response = make_response(jsonify(toret), 500)
            response.headers["Content-Type"] = "application/json"
            return response
        elif job.get_status()=='finished':
            pipeline_result = job.return_value
            if not pipeline_result[0]:
                logger.error('The job '+str(job.id)+' failed: '+str(pipeline_result[1]))
                toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta, 'pipeline_status': pipeline_result[0], 'pipeline_error_msg': pipeline_result[2]}
                logger.debug(toret)
                response = make_response(jsonify(toret), 500)
                response.headers["Content-Type"] = "application/json"
                return response
            else:
                logger.debug('Successfull run')
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
            logger.debug(toret)
            response = make_response(jsonify(toret), 102)
            response.headers["Content-Type"] = "application/json"
            return response


class RestQueryJob(Resource):
    def post(self):
        params = json.load(request.files['data'])
        redis_conn = Redis()
        logger.debug('Received a query for job: '+str(params['job_id']))
        job = Job.fetch(params['job_id'], connection=redis_conn)
        logger.debug(job)
        toret = {'job_id': job.id, 'status': job.get_status(), 'meta': job.meta}
        logger.debug(toret)
        response = make_response(jsonify(toret), 200)
        response.headers["Content-Type"] = "application/json"
        return response
        

api.add_resource(RestApp, '/REST')
api.add_resource(RestSubmitJob, '/REST/SubmitJob')
api.add_resource(RestQueryJob, '/REST/QueryJob')
api.add_resource(RestRetreiveJob, '/REST/RetreiveJob')

if __name__== "__main__":
    #handler = RotatingFileHandler('metaxime.log', maxBytes=10000, backupCount=1)
    #handler.setLevel(logging.DEBUG)
    #logger.addHandler(handler)
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=False)
