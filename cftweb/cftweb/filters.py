"""
Custom Jinja filters for use in CFT templates.

Learn more about Jinja filters at,
http://jinja.pocoo.org/docs/dev/api/#writing-filters

Also see typical usage in `templates/individuals.html` and `templates/index.html`

"""

from __future__ import print_function
import collections
from cftweb import app

# I'm going to start using cottonmouth for hiccup templating (look it up).
# The basic idea is to use basic python data structures as a way to describe html dom.
# This ends up being nice, because of instead of having to mess with the semantics of constructing valid html
# strings, we can focus on constructing simple python data structures, leaving details of closing braces etc
# for the translation lib.
from cottonmouth import html, tags
#from cottonmouth.html import render
#from cottonmouth.tags import html, head, body, title, meta, link, h1


@app.template_filter()
def clustersort(value, by=['run', 'seed', 'v_gene', 'd_gene', 'j_gene'],reverse=False):
    """Sort a dict of cluster objects by attributes supplied in `by`.

        {% for item in mydict|dictsort %}
            sort the dict by key, case insensitive

        {% for item in mydict|dictsort(true) %}
            sort the dict by key, case sensitive

        {% for item in mydict|dictsort(false, 'value') %}
            sort the dict by key, case insensitive, sorted
            normally and ordered by value.
    """
    def sort_func(item):
        value = item[1]
        return tuple(getattr(value, k) for k in by)

    return sorted(value.items(), key=sort_func, reverse=reverse)


# unique Jinja filter cribbed from ansible,
# https://github.com/ansible/ansible/blob/6787fc70a643fb6e2bdd2c6a6202072d21db72ef/lib/ansible/plugins/filter/mathstuff.py#L28
#
@app.template_filter()
def unique(a):
    if isinstance(a,collections.Hashable):
        c = set(a)
    else:
        c = []
        for x in a:
            if x not in c:
                c.append(x)
    return c


# Now for sequence annotation

# You'll see here my first usage of cottonmouth for hiccup (see notes above)
# Note also that we've pulled things apart somewhat

def mutations(seq, cluster, seq_mode="dna"):
    "Returns a set of mutation indices"
    ## This is currently just being mocked to #{}
    #seq = seq.upper()
    # random indices to mock up to 10 mutations in each sequence
    #import random
    #lenseq = len(seq)
    #if lenseq >= 10:
    #    mutations = set(random.sample(range(lenseq), random.randrange(10)))
    #else:
    #    mutations = []
    mutations = set([]) # <-- eventually read from metadata
    return mutations

def chain(seq, cluster):
    if cluster.d_start == cluster.d_end:
        chain = 'light'
    else:
        chain = 'heavy'
    return chain


# Note the old string munging technique

def annotate_dna(seq, cluster):
    foo = ''#seq[:vIndex]
    # NOTE: this will make people laugh at you, recode to hierarchical spans
    for i in range(int(cluster.v_start), (cluster.j_end)):
        if int(cluster.v_start) <= i < int(cluster.v_end):
            gene_class = 'Vgene'
        elif (chain(seq, cluster) == 'heavy' and (int(cluster.v_end) <= i < int(cluster.d_start) or int(cluster.d_end) <= i < int(cluster.j_start))) or \
             (chain(seq, cluster) == 'light' and int(cluster.v_end) <= i < int(cluster.j_start)):
            gene_class = 'N'
        elif int(cluster.d_start) <= i < int(cluster.d_end):
            gene_class = 'Dgene'
        elif int(cluster.j_start) <= i:
            gene_class = 'Jgene'
        foo += '<span class="' + gene_class + \
               (' mutation' if i in mutations(seq, cluster, 'dna') else '') + \
               (' CDR3' if int(cluster.cdr3_start) <= i < int(cluster.cdr3_start) + int(cluster.cdr3_length) else '') + \
               '"">' + seq[i] + '</span>'
    return foo


# Compared to the simplicity of constructing a data structure representing the corresponding dom?

def annotate_aa(seq, cluster):
    muts = mutations(seq, cluster, 'aa')
    hiccup = ["span.sequence.aa",
              (['span.aa-' + bp.lower() + ('.mutation' if i in muts else ''),
                    bp.upper()]
                  for i, bp in enumerate(seq))]
    return html.render(hiccup)

# Basic idea is that the first element of the list is the html tag (possibly followed by `.someclass#theid`
# for classes then ids), then optionally a dict for html attrs, then followed by child elements (where
# generators are "inlined" as children).
# Very simple mapping from data to dom => very big win.


@app.template_filter()
def annotate(seq, cluster, seq_mode="dna"):
    # Make sure seq_mode is supported
    assert(seq_mode in set(["dna", "aa"]))
    # Return if empyt seq or seq mode is aa (for now; might annotate aa eventually...)
    if len(seq) == 0: return seq

    if seq_mode == "dna":
        return annotate_dna(seq, cluster)
    else:
        return annotate_aa(seq, cluster)


