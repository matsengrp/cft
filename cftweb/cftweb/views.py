#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Description

Usage:
    fill in typical command-line usage

'''
from __future__ import print_function

import flask
from cStringIO import StringIO
from flask import render_template, Response, url_for, request
from flask_breadcrumbs import register_breadcrumb

from Bio import SeqIO

from cftweb import app
import filters

# annoying error indicator hack; semantically irrelevant; need to load filters for filters to be loaded somewhere in html templates
_ = filters

def base_renderdict(dataset_id, updates={}):
    # God damn mutability...
    renderdict = app.config['CLUSTERS'].build_info(dataset_id)
    renderdict.update(app.config['CFTWEB_BUILD_INFO'])
    renderdict['commit_url'] = "https://github.com/matsengrp/cft/tree/" + renderdict["commit"]
    renderdict['short_commit'] = renderdict["commit"][0:10]
    renderdict['app_commit_url'] = "https://github.com/matsengrp/cft/tree/" + renderdict["app_commit"]
    renderdict['app_short_commit'] = renderdict["app_commit"][0:10]
    renderdict['dataset_id'] = dataset_id
    renderdict.update(updates)
    return renderdict


# This is the more or less what the routing scheme looks like.
# At the very top level we have our dataset ids; everything nests nicely within this.
# The only exception to this is the cluster page route, which doesn't need this nesting, as the cluster id is
# computed to be unique between datasets.

#   /v8-dnaml-2017.10.23/subjects
#   /v8-dnaml-2017.10.23/seeds/QA255
#   /v8-dnaml-2017.10.23/clusters
#   /v8-dnaml-2017.10.23/clusters/QA255/006-Vh
#   /clusters/c03943943


def by_cluster_response(query):
    query = {k: v for k, v in query.iteritems() if v != None}
    clusters = app.config['CLUSTERS'].query(query)
    renderdict = base_renderdict({'clusters': clusters})
    return render_template('by_cluster.html', **renderdict)


@app.route('/')
@app.route('<dataset_id>/clusters')
@register_breadcrumb(app, '.', 'Cluster Index')
def index(dataset_id=None):
    # For now, just pick some random build; Would be nice to troll through the build info and pick something
    # more sensible based on date and dnaml/dnapars preference.
    dataset_id = dataset_id or app.config['CLUSTERS'].build_info.keys()[0]
    return by_cluster_response({'dataset_id': dataset_id})


@app.route("/<dataset_id>/subjects")
@register_breadcrumb(app, '.', 'By Subject')
def by_subject(**params):
    clusters = app.config['CLUSTERS'].query(params)
    renderdict = base_renderdict({'clusters': clusters})
    return render_template('by_subject.html', **renderdict)


def subject_id_bcrumb():
    subject_id = request.view_args['subject_id']
    return [{'text': subject_id, 'url': url_for('by_timepoint', subject_id=subject_id)}]

@app.route("/<dataset_id>/timepoints/<subject_id>")
@register_breadcrumb(app, '.subject_id', '<ignored>',
                     dynamic_list_constructor=subject_id_bcrumb)
def by_timepoint(**params):
    params['clusters'] = app.config['CLUSTERS'].query(params)
    renderdict = base_renderdict(params)
    return render_template('by_timepoint.html', **renderdict)


def timepoint_bcrumb():
    subject_id = request.view_args['subject_id']
    timepoint = request.view_args['timepoint']
    return [{'text': timepoint, 'url': url_for('by_seed', subject_id=subject_id, timepoint=timepoint)}]

@app.route("/<dataset_id>/seeds/<subject_id>/<timepoint>")
@register_breadcrumb(app, '.subject_id.timepoint', '<ignored>',
                     dynamic_list_constructor=timepoint_bcrumb)
def by_seed(**params):
    params['clusters'] = app.config['CLUSTERS'].query(params)
    renderdict = base_renderdict(params)
    return render_template('by_seed.html', **renderdict)


def seed_bcrumb():
    subject_id = request.view_args['subject_id']
    timepoint = request.view_args['timepoint']
    seedid = request.view_args['seedid']
    return [{'text': seedid, 'url': url_for('by_seed', subject_id=subject_id, timepoint=timepoint, seedid=seedid)}]

@app.route("/<dataset_id>/clusters/<subject_id>/<timepoint>/<seedid>")
@register_breadcrumb(app, '.subject_id.timepoint.seed', '<ignored>',
                     dynamic_list_constructor=seed_bcrumb)
def by_cluster(**params):
    return by_cluster_response(params)


@app.route("/cluster/<id>")
def cluster_page(id=None):
    cluster = app.config['CLUSTERS'].get_by_id(id)
    clustering_step_siblings = app.config['CLUSTERS'].clustering_step_siblings(id)

    # We are now using the lineage kw_args param to get the lineage we want to display
    seed_name = cluster.seed
    focus_nodes = request.args.get('focus_nodes', seed_name).split(",")

    seq_mode = request.args.get('seq_mode', 'aa') # also accepts 'dna' for dna; defaults to amino acid

    renderdict = base_renderdict({
        'cluster': cluster,
        'clustering_step_siblings': clustering_step_siblings,
        'focus_nodes': focus_nodes,
        'lineage_seqs': cluster.multi_lineage_seqs(focus_nodes, seq_mode),
        #'lineage_seqs': cluster.lineage_seqs(focus_nodes, seq_mode),
        'seed_seq': cluster.seed_seq(seq_mode),
        'seq_mode': seq_mode,
        'svg': cluster.svgstr()})

    return render_template('cluster.html', **renderdict)

def to_fasta(seqs):
    fp = StringIO()
    SeqIO.write(seqs, fp, 'fasta')
    return fp.getvalue()


@app.route("/cluster/<id>/lineage/<focus_nodes>.lineage.<seq_mode>.fa")
def lineage_download(id, focus_nodes, seq_mode):
    cluster = app.config['CLUSTERS'].get_by_id(id)
    seqs = cluster.lineage_seqs(focus_nodes, seq_mode)

    fasta = to_fasta(seqs)
    return Response(fasta, mimetype="application/octet-stream")


@app.route("/cluster/sequences/<id>.fa")
def cluster_fasta(id=None):
    cluster = app.config['CLUSTERS'].get_by_id(id)

    fasta = to_fasta(cluster.sequences())
    return Response(fasta, mimetype="application/octet-stream")


# TODO We may want to make seedlineage a function
@app.route("/download/sequences/<id>.seedlineage.fa")
def seedlineage_fasta(id=None):
    cluster = app.config['CLUSTERS'].get_by_id(id)
    seqs = cluster.lineage_seqs(cluster.seed)

    fasta = to_fasta(seqs)
    return Response(fasta, mimetype="application/octet-stream")


@app.route("/download/sequences/<id>_aa.fa")
def cluster_aa_fasta(id=None):
    cluster = app.config['CLUSTERS'].get_by_id(id)
    return flask.send_file(cluster.cluster_aa)

@app.route("/download/sequences/<id>.seedlineage_aa.fa")
def seedlineage_aa_fasta(id=None):
    cluster = app.config['CLUSTERS'].get_by_id(id)
    seqs = cluster.lineage_seqs(cluster.seed, seq_mode='aa')

    fasta = to_fasta(seqs)
    return Response(fasta, mimetype="application/octet-stream")


@app.route("/cluster/tree/<id>.svg")
def svg_view(id=None):
    # TODO should just have middleware for passing through the cluster
    cluster = app.config['CLUSTERS'].get_by_id(id)
    return flask.send_file(cluster.svg)


