import os

import itertools
import logging
import json
from flask import Flask, g
from flask_breadcrumbs import Breadcrumbs

from cluster import Cluster

import os.path
os.environ['CFTWEB_SETTINGS'] = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), 'config.py')

app = Flask(__name__)
app.config.from_envvar('CFTWEB_SETTINGS')
app.config['DEBUG'] = True

# Initialize Flask-Breadcrumbs
Breadcrumbs(app=app)

# iterate through json files under`dir`
def content_file_iterator(dir):
    for root, dirs, files in os.walk(dir):
        for f in files:
            if os.path.splitext(f)[1] == '.json':
                yield os.path.join(root, f)


# Initiate engine before the first request
@app.before_first_request
def before_first_request():
    options = app.config['OPTIONS']

    if options.dir:
        objects = [Cluster.fromfile(f) for f in content_file_iterator(options.dir)]
    else:
        objects = [Cluster.fromfile(options.file)]
        
    objects = itertools.chain.from_iterable(objects)
    objects = [o for o in objects if o is not None]
    objects = dict([(g.id, g) for g in objects])
    app.config['CLUSTERS'] = objects

    # Configure logging
    if not app.debug:
        from logging.handlers import SMTPHandler, RotatingFileHandler
        from logging import Formatter
        mail_handler = SMTPHandler('127.0.0.1', 'igdbweb@fredhutch.org',
                                   app.config['ADMINS'], 'igdb failed')
        mail_handler.setLevel(logging.ERROR)
        log_path = app.config['LOG_PATH']
        file_handler = RotatingFileHandler(
            log_path, maxBytes=1 << 20, backupCount=5)
        file_handler.setFormatter(
            Formatter('%(asctime)s %(levelname)s: %(message)s '
                      '[in %(pathname)s:%(lineno)d]'))
        app.logger.setLevel(logging.INFO)
        file_handler.setLevel(logging.INFO)
        app.logger.addHandler(file_handler)
        app.logger.addHandler(mail_handler)


# @app.errorhandler(Exception)
# def all_exception_handler(error):
#    return 'Error', 500

import views
