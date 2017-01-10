#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Convert CFT json files into a static website for data exploration.
'''
from __future__ import print_function

import argparse

from cftweb import app

# pass commandline arguments to flask app
# http://flask.pocoo.org/snippets/133/

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='python -m cftbweb', description=__doc__)

    parser.add_argument(
        '-f',
        '--file',
        help="""Name of metadata file: [default "%(default)s"]""")

    parser.add_argument(
        '-H',
        '--host',
        default="0.0.0.0",
        help='Hostname of the Flask app [%(default)s]')
    parser.add_argument(
        '-P',
        '--port',
        default=5000,
        type=int,
        help='Port for the Flask app  [%(default)d]')
    parser.add_argument(
        '-d',
        '--debug',
        default=False,
        action="store_true",
        help='turn on debugging output (and turn off email, slack and file logging)')
    parser.add_argument(
        '-e',
        '--email',
        default=False,
        action="store_true",
        help="send error message notifications to admins (see config) via email")
    parser.add_argument(
        '-s',
        '--slack',
        default=False,
        action="store_true",
        help="send error message notifications to slack")
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=False,
        help='verbose output')
    

    global options
    options = parser.parse_args()
    app.config['OPTIONS'] = options

    if not options.debug:
        try:
            from gevent import wsgi
            print("Running with WSGI server")
            http_server = wsgi.WSGIServer((options.host, options.port), app)
            http_server.serve_forever()
        except ImportError:
            print("Unable to run with WSGI server; using built in server (may crash spontaneously)")
            print("See: http://stackoverflow.com/questions/37962925/flask-app-get-ioerror-errno-32-broken-pipe")
            print("Run in debug mode to prevent this warning")
            app.run(debug=options.debug, host=options.host, port=options.port)

    else:
        app.run(debug=options.debug, host=options.host, port=int(options.port))

    
