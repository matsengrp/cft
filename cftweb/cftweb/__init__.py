import os

import logging
import json
import sys
import datetime
import getpass
import subprocess
from flask import Flask
from flask_breadcrumbs import Breadcrumbs

from cluster import Cluster

import os.path
os.environ['CFTWEB_SETTINGS'] = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), 'config.py')

app = Flask(__name__)
app.config.from_envvar('CFTWEB_SETTINGS')
app.config['DEBUG'] = True

def git(*args):
    return subprocess.check_output(['git'] + list(args))

app.config['CFTWEB_BUILD_INFO'] = {
        'app_date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'app_command': " ".join(sys.argv),
        'app_workdir': os.getcwd(),
        'app_user': getpass.getuser(),
        'app_commit': git('rev-parse', 'HEAD'),
        'app_status': git('status', '--porcelain', '.')}

# Initialize Flask-Breadcrumbs
Breadcrumbs(app=app)


def email_on_error(app):
    from logging.handlers import SMTPHandler
    mail_handler = SMTPHandler('127.0.0.1', 'csmall@fredhutch.org',
                               app.config['ADMINS'], 'cftweb failed')
    mail_handler.setLevel(logging.ERROR)
    app.logger.addHandler(mail_handler)


def slack_on_error(app):
    import logging
    from slacker_log_handler import SlackerLogHandler
    slack_handler = SlackerLogHandler(os.environ['SLACK_TOKEN'],
                                      'cft',
                                      stack_trace=True,
                                      username='cftweb-robot')
    # Should probably set everything to use the app.config['LOGGING_FORMAT'] as shown here
    #slack_handler.setFormatter(Formatter(app.config['LOGGING_FORMAT']))
    slack_handler.setLevel(logging.ERROR)
    app.logger.addHandler(slack_handler)


def log_to_file(app):
    from logging.handlers import RotatingFileHandler
    from logging import Formatter
    log_path = app.config['LOG_PATH']
    file_handler = RotatingFileHandler(
        log_path, maxBytes=1 << 20, backupCount=5)
    file_handler.setFormatter(
        Formatter('%(asctime)s %(levelname)s: %(message)s '
                  '[in %(pathname)s:%(lineno)d]'))
    file_handler.setLevel(logging.INFO)
    app.logger.addHandler(file_handler)



# Initiate engine before the first request
@app.before_first_request
def before_first_request():
    options = app.config['OPTIONS']

    objects = Cluster.fromfile(options.file)

    objects = [o for o in objects if o is not None]
    objects = dict((g.id, g) for g in objects)
    app.config['CLUSTERS'] = objects
    with open(options.file) as fh:
        app.config['DATA_BUILD_INFO'] = json.load(fh)["build_info"]

    # Configure logging
    if not app.debug:
        # Turning this off for now...
        log_to_file(app)
        if app.config['OPTIONS'].email:
            email_on_error(app)
        if app.config['OPTIONS'].slack:
            slack_on_error(app)
        app.logger.setLevel(logging.INFO)


# @app.errorhandler(Exception)
# def all_exception_handler(error):
#    return 'Error', 500

import views
