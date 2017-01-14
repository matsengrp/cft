#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Description

Usage:
    fill in typical command-line usage

'''
from __future__ import print_function

import os
import copy
import flask
from cStringIO import StringIO
from flask import render_template, Response, url_for, request
from flask_breadcrumbs import register_breadcrumb

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from cftweb import app
import filters

# annoying error indicator hack; semantically irrelevant; need to load filters for filters to be loaded somewhere in html templates
_ = filters

def base_renderdict(updates={}):
    # God damn mutability...
    renderdict = copy.copy(app.config['DATA_BUILD_INFO'])
    renderdict.update(app.config['CFTWEB_BUILD_INFO'])
    renderdict['commit_url'] = "https://github.com/matsengrp/cft/tree/" + renderdict["commit"]
    renderdict['short_commit'] = renderdict["commit"][0:10]
    renderdict['app_commit_url'] = "https://github.com/matsengrp/cft/tree/" + renderdict["app_commit"]
    renderdict['app_short_commit'] = renderdict["app_commit"][0:10]
    renderdict.update(updates)
    return renderdict


@app.route('/index')
@app.route("/")
@register_breadcrumb(app, '.', 'Index')
def index():
    clusters = app.config['CLUSTERS'].values()
    clusters = sorted(clusters, key=lambda c: (c.pid, c.timepoint, c.seed, c.clustering_step))
    renderdict = base_renderdict({'clusters': clusters, 'heading': "All clusters"})
    return render_template('by_cluster.html', **renderdict)


@app.route("/individuals.html")
@register_breadcrumb(app, '.', 'By Patient')
def by_pid():
    clusters = app.config['CLUSTERS'].values()
    renderdict = base_renderdict({'clusters': clusters})
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
    renderdict = base_renderdict({
        'pid' : pid,
        'clusters': clusters})
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
    renderdict = base_renderdict({
        'pid' : pid,
        'timepoint' : timepoint,
        'clusters': clusters})
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
    clusters = sorted(clusters, key=lambda c: (c.seed, c.clustering_step))

    renderdict = base_renderdict({
        'pid' : pid,
        'timepoint' : timepoint,
        'seedid' : seedid,
        'clusters': clusters})

    return render_template('by_cluster.html', **renderdict)


def lineage_seqs(cluster, focus_node_name, seq_mode="dna"):
    assert(seq_mode in set(["dna", "aa"]))
    tree = cluster.tree()
    focus_node = tree.search_nodes(name=focus_node_name)[0]
    lineage = [focus_node] + focus_node.get_ancestors()
    fname = cluster.fasta if seq_mode == "dna" else cluster.cluster_aa
    with open(fname, "rU") as fh:
        cluster_seqs = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
        #import pdb; pdb.set_trace()
        seqs = [cluster_seqs[n.name] for n in lineage if n.name in cluster_seqs]
        # I think to add in naive we have to do something like
        #naive_node = find_node(tree, '.*naive.*')
        # Then get the node name from that then do the below (with correct node name)
        #seqs.append(cluster_seqs['naive'])
    return seqs


@app.route("/cluster/<id>/cluster.html")
def cluster_page(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    # We are now using the lineage kw_args param to get the lineage we want to display
    seed_name = cluster.seed
    focus_node = request.args.get('focus_node', seed_name)

    seq_mode = request.args.get('seq_mode', 'aa') # also accepts 'dna' for dna; defaults to amino acid

    renderdict = base_renderdict({
        'cluster': cluster,
        'focus_node': focus_node,
        'lineage_seqs': lineage_seqs(cluster, focus_node, seq_mode),
        'seq_mode': seq_mode,
        'svg': cluster.svgstr()})

    return render_template('cluster.html', **renderdict)


@app.route("/download/fasta/<id>.fa")
def cluster_fasta(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    # It would probably be better just to transmit the self.fasta attr...
    def to_fasta(seqs):
        fp = StringIO()
        SeqIO.write(seqs, fp, 'fasta')
        return fp.getvalue()

    fasta = to_fasta(cluster.sequences())
    return Response(fasta, mimetype="application/octet-stream")


@app.route("/download/fasta/<id>.seedlineage.fa")
def seedlineage_fasta(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    return flask.send_file(cluster.seedlineage)


@app.route("/download/fasta/<id>_aa.fa")
def cluster_aa_fasta(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    return flask.send_file(cluster.cluster_aa)

@app.route("/download/fasta/<id>.seedlineage_aa.fa")
def seedlineage_aa_fasta(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    return flask.send_file(cluster.seedlineage_aa)



@app.route("/cluster/<id>/tree.svg")
def svg_view(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    return flask.send_file(cluster.svg)


# This is brokerz... should remove if no one needs it...
@app.route("/download/phylip/<id>.phy")
def cluster_phylip(id=None):
    clusters = app.config['CLUSTERS']
    cluster = clusters[id]

    def to_phylip(seqs):
        fp = StringIO()
        for seq in seqs:
            seq = seq.replace('.', '-')
            record = SeqRecord(seq)
            SeqIO.write(record, fp, 'phylip-relaxed')
        return fp.getvalue()

    phylip = to_phylip(cluster.sequences())

    return Response(phylip, mimetype="application/octet-stream")

