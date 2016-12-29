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

    app.run(debug=options.debug, host=options.host, port=int(options.port))


    
