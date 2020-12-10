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
from datetime import datetime
from flask import Flask, request, jsonify, send_file, abort, Response, make_response
from flask_restful import Resource, Api
from logging.handlers import RotatingFileHandler

from logging.config import dictConfig

from rq import Connection, Queue
from redis import Redis

import pipeline

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
    """Class containing the REST requests for the pipeline
    """
    def post(self):
        """Make the REST request using the POST method

        :rtype: Response
        :return: Flask Response object
        """
        params = json.load(request.files['data'])
        app.logger.debug('name: '+str(params['name']))
        app.logger.debug('SMILES: '+str(params['smiles']))
        app.logger.debug('GEM: '+str(params['gem']))
        app.logger.debug('Steps: '+str(params['steps']))
        ##### REDIS ##############
        conn = Redis()
        q = Queue('default', connection=conn, default_time_out='24h')
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            fi_name = params['name']+'__'+params['gem']+'__'+str(time.time())+'.rpcol'
            app.logger.debug('fi_name: '+str(fi_name))
            rpcoll = os.path.join(tmp_output_folder, fi_name)
            app.logger.debug('rpcoll: '+str(rpcoll))
            app.logger.debug('pipeline: '+str(pipeline.pipeline))
            async_results = q.enqueue(pipeline.pipeline, rpcoll, params['smiles'], params['gem'], int(params['steps']))
            result = None
            while result is None:
                result = async_results.return_value
                app.logger.debug('result: '+str(result))
                if async_results.get_status()=='failed':
                    return Response('Redis job failed \n '+str(result), status=500)
                time.sleep(2.0)
            app.logger.debug(result)
            if status:
                response = make_response(send_file(rpcoll, as_attachment=True, attachment_filename=fi_name, mimetype='application/x-tar'))
                response.headers['status_message'] = err_type
                return response
            else:
                app.logger.error('There was an error running the pipeline: '+str(err_type))
                return Response('There was an error running the pipeline: '+str(err_type), status=500)

api.add_resource(RestApp, '/REST')
api.add_resource(RestQuery, '/REST/Query')

if __name__== "__main__":
    handler = RotatingFileHandler('metaxime.log', maxBytes=10000, backupCount=1)
    handler.setLevel(logging.DEBUG)
    app.logger.addHandler(handler)
    app.run(host="0.0.0.0", port=8888, debug=True, threaded=True)
