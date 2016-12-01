
import re
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


# in-order iterator nodes in ETE tree.
# left subtree, then parent, then right subtree.
# Use this to order sequence names in the same order as nodes in a tree.
def iter_inorder(t):
    c = t.children
    if c:
        for x in iter_inorder(c[0]):
            yield x
        yield t
        if len(c) > 1:
            for x in iter_inorder(c[1]):
                yield x
    else:
        yield t


def iter_names_inorder(t):
    for n in iter_inorder(t):
        yield n.name

def find_node(tree, pattern):
    regex = re.compile(pattern).search
    nodes =  [ node for node in tree.traverse() for m in [regex(node.name)] if m]
    if nodes:
        # if len(nodes) > 1:
        #     warn("multiple nodes found; using first one.\nfound: {}".format([n.name for n in nodes]))
        return nodes[0]

def fake_seq():
    """
    Return an empty placeholder sequence with an empty name and an empty
    sequence.  This is used ot get the spacing right when lining up
    sequences adjacent to ascii-art tree topologies
    """
    record = SeqRecord(Seq("",
                   IUPAC.unambiguous_dna),
                   id="", name="",
                   description="placeholder blank sequence")
    return record


def sort_tree(tree, direction=0, predicate=lambda x: True):
    """
    Sort the nodes in a tree according to the supplied predicate function.

    This is a generic solution to the problem is showing the 'seed' node
    at the top of the tree.

    Unlike a comparator function used in most sort routines, the
    predicate function used here returns only `true` or `false`.  It
    should return `true` if the node should be to the left of other
    nodes, and false otherwise.

    Typical usage is
        sort_tree(tree, direction=1, predicate=lambda n: 'seed' in n.name)
    to move all nodes with the word 'seed' in their name to the right
    (because direction == 1).  The predicate function could also use a
    regular expression to match the node name, or could depend on a
    different attribute altogether.  The predicate function may return
    true for multiple nodes, the indicated nodes will move to one side
    of their respective subtree, but the order of nodes within the
    indicated group is undefined.

    if `direction` is 1, the direction of movement is reversed.
    """

    pred = False
    if not tree.is_leaf():
        n2s = {}
        for n in tree.get_children():
            s = sort_tree(n, direction=direction, predicate=predicate)
            n2s[n] = s

        tree.children.sort(key=lambda x: n2s[x], reverse=(direction==1))
        pred = any(n2s.values())
    else:
        pred = predicate(tree)

    return pred
