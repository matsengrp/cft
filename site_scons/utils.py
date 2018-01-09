import copy

# Utility functions
# -----------------

def merge_dicts(d1, d2):
    d = copy.deepcopy(d1)
    d.update(d2)
    return d


def get_in(d, ks, default=None):
    if len(ks) > 1:
        return get_in(d.get(ks[0], {}), ks[1:], default=default)
    else:
        return d.get(ks[0], default)
