

from __future__ import print_function


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

def register(app):
    """
    Register all filters with an application
    """
    app.jinja_env.filters['clustersort'] = clustersort
