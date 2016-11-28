
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
    """Return a random placeholder sequence with an empty name and an empty sequence.
    This is used ot get the spacing right when lining up sequences adjacent to ascii-art tree topologies
    """
    record = SeqRecord(Seq("",
                   IUPAC.unambiguous_dna),
                   id="", name="",
                   description="placeholder blank sequence")
    return record

