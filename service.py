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

import logging
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort, Response, make_response
from flask_restful import Resource, Api
from logging.handlers import RotatingFileHandler


from rq import Connection, Queue
from redis import Redis

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
class RestQuery(Resource):
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
        conn = Redis()
        q = Queue('default', connection=conn, default_time_out='24h')
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            fi_name = params['name']+'__'+params['gem']+'__'+str(time.time())+'.rpcol'
            logger.debug('fi_name: '+str(fi_name))
            rpcoll = os.path.join(tmp_output_folder, fi_name)
            logger.debug('rpcoll: '+str(rpcoll))
            logger.debug('pipeline: '+str(pipeline.pipeline))
            async_results = q.enqueue(pipeline.pipeline, str(rpcoll), str(params['smiles']), str(params['gem']), int(params['steps']))
            result = None
            while result is None:
                result = async_results.return_value
                if async_results.get_status()=='failed':
                    return Response('Redis job failed \n '+str(result), status=500)
                time.sleep(5.0)
            if status:
                response = make_response(send_file(rpcoll, as_attachment=True, attachment_filename=fi_name, mimetype='application/x-tar'))
                response.headers['status_message'] = err_type
                return response
            else:
                logger.error('There was an error running the pipeline: '+str(err_type))
                return Response('There was an error running the pipeline: '+str(err_type), status=500)

api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')

if __name__== "__main__":
    #handler = RotatingFileHandler('metaxime.log', maxBytes=10000, backupCount=1)
    #handler.setLevel(logging.DEBUG)
    #logger.addHandler(handler)
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=True)
