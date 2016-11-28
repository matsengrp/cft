#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Description

Usage:
    fill in typical command-line usage

'''
from __future__ import print_function

import os
import sys

import cluster
import datetime
import getpass
import hashlib
from cStringIO import StringIO
from jinja2 import Environment, FileSystemLoader
from flask import render_template, abort, Response, url_for, request, g, jsonify
from flask_breadcrumbs import register_breadcrumb

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from cftweb import app
import filters

@app.route('/index')
@app.route("/")
@register_breadcrumb(app, '.', 'Index')
def index():
    clusters = app.config['CLUSTERS'].values()
    renderdict = {
        'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'command': " ".join(sys.argv),
        'workdir': os.getcwd(),
        'user': getpass.getuser(),
        'title': 'CFT Cluster Visualization',
        'clusters': clusters
    }

    return render_template('by_cluster.html', **renderdict)

@app.route("/individuals.html")
@register_breadcrumb(app, '.', 'By Patient')
def by_pid():
    clusters = app.config['CLUSTERS'].values()
    renderdict = {
        'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'command': " ".join(sys.argv),
        'workdir': os.getcwd(),
        'user': getpass.getuser(),
        'title': 'CFT Cluster Visualization',
        'clusters': clusters
    }

    return render_template('by_pid.html', **renderdict)

def pid_bc():
    pid = request.view_args['pid']
    return [{'text': pid, 'url': url_for('by_timepoint', pid=pid)}]


@app.route("/<pid>/timepoints.html")
@register_breadcrumb(app, '.pid', '<ignored>',
                     dynamic_list_constructor=pid_bc)
def by_timepoint(pid):
    clusters = app.config['CLUSTERS'].values()
    clusters = [c for c in clusters if c.pid == pid]
    renderdict = {
        'pid' : pid,
        'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'command': " ".join(sys.argv),
        'workdir': os.getcwd(),
        'user': getpass.getuser(),
        'title': 'CFT Cluster Visualization',
        'clusters': clusters
    }

    return render_template('by_timepoint.html', **renderdict)

def timepoint_bc():
    pid = request.view_args['pid']
    timepoint = request.view_args['timepoint']
    return [{'text': timepoint, 'url': url_for('by_seed', pid=pid, timepoint=timepoint)}]

@app.route("/<pid>/<timepoint>/seeds.html")
@register_breadcrumb(app, '.pid.timepoint', '<ignored>',
                     dynamic_list_constructor=timepoint_bc)
def by_seed(pid, timepoint):
    clusters = app.config['CLUSTERS'].values()
    clusters = [c for c in clusters if c.pid == pid]
    clusters = [c for c in clusters if c.timepoint == timepoint]

    renderdict = {
        'pid' : pid,
        'timepoint' : timepoint,

        'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'command': " ".join(sys.argv),
        'workdir': os.getcwd(),
        'user': getpass.getuser(),
        'title': 'CFT Cluster Visualization',
        'clusters': clusters
    }

    return render_template('by_seed.html', **renderdict)

def seed_bc():
    pid = request.view_args['pid']
    timepoint = request.view_args['timepoint']
    seedid = request.view_args['seedid']
    return [{'text': seedid, 'url': url_for('by_seed', pid=pid, timepoint=timepoint, seedid=seedid)}]

@app.route("/<pid>/<timepoint>/<seedid>/clusters.html")
@register_breadcrumb(app, '.pid.timepoint.seed', '<ignored>',
                     dynamic_list_constructor=seed_bc)
def by_cluster(pid, timepoint, seedid):
    clusters = app.config['CLUSTERS'].values()
    clusters = [c for c in clusters if c.pid == pid]
    clusters = [c for c in clusters if c.timepoint == timepoint]
    clusters = [c for c in clusters if c.seedid == seedid]

    renderdict = {
        'pid' : pid,
        'timepoint' : timepoint,
        'seedid' : seedid,

        'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        'command': " ".join(sys.argv),
        'workdir': os.getcwd(),
        'user': getpass.getuser(),
        'title': 'CFT Cluster Visualization',
        'clusters': clusters
    }

    return render_template('by_cluster.html', **renderdict)


@app.route("/cluster/<id>/overview.html")
def cluster_svgpage(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    renderdict = {
        'id': id,
        'svg': cluster.svgstr(),
        'records': cluster.sequences(),
    }

    return render_template('cluster.html', **renderdict)

@app.route("/cluster/<id>/tree.html")
def cluster_tree(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]
    renderdict = {'svg': cluster.svgstr(), }

    return render_template('tree.html', **renderdict)


@app.route("/cluster/<id>/sequences.html")
def cluster_sequences(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    renderdict = {'records': cluster.sequences(), }

    return render_template('sequences.html', **renderdict)


@app.route("/download/fasta/<id>.fa")
def cluster_fasta(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    def to_fasta(seqs):
        fp = StringIO()
        SeqIO.write(seqs, fp, 'fasta')
        return fp.getvalue()

    fasta = to_fasta(cluster.sequences())
    return Response(fasta, mimetype="application/octet-stream")


@app.route("/download/phylip/<id>.phy")
def cluster_phylip(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    def to_phylip(seq):
        fp = StringIO()
        for seq in seqs:
            seq = seq.replace('.', '-')
            record = SeqRecord(seq)
            SeqIO.write(record, fp, 'phylip-relaxed')
        return fp.getvalue()

    phylip = to_phylip(cluster.sequences())

    return Response(phylip, mimetype="application/octet-stream")

@app.route("/cluster/<id>/asciiart.html")
def cluster_page(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    renderdict = {
        'cluster': cluster,
        'asciiart': cluster.tree().get_ascii(),
        'records': cluster.sequences(),
    }

    return render_template('asciiart.html', **renderdict)


