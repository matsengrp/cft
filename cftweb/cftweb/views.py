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

# Utility function for subsetting dicts to create minimal cluster query and routing params
def dict_subset(d, attrs):
    return {k: d[k] for k in attrs}

# Applies above utility function towards request.view_args for query/routing params and such
def get_view_args(attrs):
    cluster_id = request.view_args.get('id')
    params = app.config['CLUSTERS'].get_by_id(cluster_id).__dict__ if cluster_id else request.view_args
    return dict_subset(params, attrs)

# Default dictionary of data to render; should probably be cleaned up a bit, since the build info could be
# left nested under a single attr pretty easily (since only used in one place, and since can change between
# datasets, feels weird)
def base_renderdict(params):
    renderdict = app.config['CLUSTERS'].build_info(params['dataset_id'])
    renderdict['datasets'] = app.config['CLUSTERS']._build_info.keys()
    print("datasets:", renderdict['datasets'])
    renderdict.update(app.config['CFTWEB_BUILD_INFO'])
    renderdict['commit_url'] = "https://github.com/matsengrp/cft/tree/" + renderdict["commit"]
    renderdict['short_commit'] = renderdict["commit"][0:10]
    renderdict['app_commit_url'] = "https://github.com/matsengrp/cft/tree/" + renderdict["app_commit"]
    renderdict['app_short_commit'] = renderdict["app_commit"][0:10]
    renderdict.update(params)
    return renderdict


# ROUTING SCHEME!
# ===============

# This is the more or less what the routing scheme looks like.
# At the very top level we have our dataset ids; everything nests nicely within this.
# The only exception to this is the cluster page route, which doesn't need this nesting, as the cluster id is
# computed to be unique between datasets.

#   /v8-dnaml-2017.10.23/subjects
#   /v8-dnaml-2017.10.23/seeds/QA255
#   /v8-dnaml-2017.10.23/clusters
#   /v8-dnaml-2017.10.23/clusters/QA255/006-Vh
#   /clusters/c03943943


def by_cluster_response(params):
    params = {k: v for k, v in params.iteritems() if v != None}
    renderdict = base_renderdict(params)
    renderdict['clusters'] = app.config['CLUSTERS'].query(params)
    return render_template('by_cluster.html', **renderdict)


# Cluster index
# -------------

def clusters_bcrumb():
    params = get_view_args(['dataset_id'])
    return [{'text': params['dataset_id'] + ' clusters', 'url': url_for('index', **params)}]

@app.route('/<dataset_id>/clusters')
@register_breadcrumb(app, '.clusterindex', 'Cluster Index',
        dynamic_list_constructor=clusters_bcrumb)
def index(dataset_id=None):
    return by_cluster_response({'dataset_id': dataset_id})

@app.route('/')
def home():
    # For now, just pick some random build; Would be nice to troll through the build info and pick something
    # more sensible based on date and dnaml/dnapars preference.
    dataset_id = app.config['CLUSTERS']._build_info.keys()[0]
    return flask.redirect(url_for('index', dataset_id=dataset_id))


# Subjects index
# --------------

def subjects_bcrumb():
    params = get_view_args(['dataset_id'])
    return [{'text': params['dataset_id'] + ' subjects', 'url': url_for('by_subject', **params)}]

@app.route("/<dataset_id>/subjects")
@register_breadcrumb(app, '.subjects', 'By Subject',
        dynamic_list_constructor=subjects_bcrumb)
def by_subject(**params):
    params['clusters'] = app.config['CLUSTERS'].query(params)
    return render_template('by_subject.html', **base_renderdict(params))


# Subject > seeds
# ---------------------------

def subject_bcrumb():
    params = get_view_args(['dataset_id', 'subject_id'])
    return [{'text': params['subject_id'], 'url': url_for('by_subject', **params)}]

@app.route("/<dataset_id>/seeds/<subject_id>")
@register_breadcrumb(app, '.subjects.seeds', '',
                     dynamic_list_constructor=subject_bcrumb)
def by_seed(**params):
    params['clusters'] = app.config['CLUSTERS'].query(params)
    renderdict = base_renderdict(params)
    return render_template('by_seed.html', **renderdict)


# Subject > seed > cluster steps
# ------------------------------------------

def seed_bcrumb():
    params = get_view_args(['dataset_id', 'subject_id', 'seedid'])
    return [{'text': params['seedid'], 'url': url_for('by_seed', **params)}]

@app.route("/<dataset_id>/clusters/<subject_id>/<seedid>")
@register_breadcrumb(app, '.subjects.seeds.clusters', '',
                     dynamic_list_constructor=seed_bcrumb)
def by_cluster(**params):
    return by_cluster_response(params)


# Cluster page
# ------------

def cluster_bcrumb():
    cluster_id = request.view_args['id']
    cluster = app.config['CLUSTERS'].get_by_id(cluster_id)
    return [{'text': cluster.clustering_step, 'url': url_for('cluster_page', id=cluster_id)}]

@app.route("/cluster/<id>")
@register_breadcrumb(app, '.subjects.seeds.clusters.cluster', '',
                     dynamic_list_constructor=cluster_bcrumb)
def cluster_page(id=None):
    cluster = app.config['CLUSTERS'].get_by_id(id)
    clustering_step_siblings = app.config['CLUSTERS'].clustering_step_siblings(id)

    # We are now using the lineage kw_args param to get the lineage we want to display
    seed_name = cluster.seed
    focus_nodes = request.args.get('focus_nodes', seed_name).split(",")

    seq_mode = request.args.get('seq_mode', 'aa') # also accepts 'dna' for dna; defaults to amino acid

    renderdict = base_renderdict({
        'cluster': cluster,
        'dataset_id': cluster.dataset_id,
        'clustering_step_siblings': clustering_step_siblings,
        'focus_nodes': focus_nodes,
        'lineage_seqs': cluster.multi_lineage_seqs(focus_nodes, seq_mode),
        #'lineage_seqs': cluster.lineage_seqs(focus_nodes, seq_mode),
        'seed_seq': cluster.seed_seq(seq_mode),
        'seq_mode': seq_mode,
        'svg': cluster.svgstr()})

    return render_template('cluster.html', **renderdict)


# Secondary cluster routes/handlers (downloads etc)
# -------------------------------------------------

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


