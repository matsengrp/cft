#! /bin/env python

import argparse

def outfile2seqs(outfile='outfile'):
    """
    give me a phylip dnaml outfile and I''ll give you a dictionary of sequences,
    including ancestral, and a dictionary of parent assignments with branch length
    """
    # parse all sequences from phylip outfile
    outfiledat = [block.split('\n\n\n')[0].split('\n\n')[1:] for block in open(outfile, 'r').read().split('Probable sequences at interior nodes:\n')[1:]]

    # only one tree in outfile
    assert len(outfiledat) == 1

    sequences = {}
    for j, block in enumerate(outfiledat[0]):
        if j == 0:
            for line in block.split('\n'):
                fields = line.split()
                if len(fields) == 0:
                    continue
                name = fields[0]
                seq = ''.join(fields[1:])
                sequences[name] = seq
        else:
            for line in block.split('\n'):
                fields = line.split()
                if len(fields) == 0: continue
                name = fields[0]
                seq = ''.join(fields[1:])
                sequences[name] += seq

    # now a second pass to get the parent data
    outfiledat = [block.split('\n\n\n')[0].split('\n\n')[1] for block in open(outfile, 'r').read().split('Between        And            Length      Approx. Confidence Limits\n')[1:]]

    # only one tree in outfile
    assert len(outfiledat) == 1

    parents = {}
    for line in outfiledat[0].split('\n'):
        fields = line.split()
        if len(fields) == 0:
            continue
        parent, child = fields[0:2]
        distance = float(fields[2])
        assert child not in parents
        parents[child] = (parent, distance)

    return sequences, parents

def main():
    from ete3 import PhyloTree
    parser = argparse.ArgumentParser(description='give me a phylip dnaml outfile and I''ll give you an alignment (including ancestral sequences) and a tree')
    parser.add_argument('--outfile', type=str, default='outfile', help='dnaml outfile (verbose output with inferred ancestral sequences, option 5). Perhaps confusingly, this is the input file to this program')
    parser.add_argument('--naive', type=str, default= None, help='name of naive (germline) sequence')
    parser.add_argument('--seed', type=str, default= None, help='name of seed sequence')
    args = parser.parse_args()

    sequences, parents = outfile2seqs(args.outfile)

    # write alignment to fasta
    aln = args.outfile+'.tree.fa'
    with open(aln, 'w') as f:
        for name in sequences:
            f.write('>'+name+'\n')
            f.write(sequences[name]+'\n')

    # build an ete tree
    # first a dictionary of disconnected nodes
    nodes = {name:PhyloTree(name=name, dist=parents[name][1] if name in parents else None, format=1) for name in sequences}
    # connect the nodes using the parent data
    orphan_nodes = 0
    for name in sequences:
        if name in parents:
            nodes[parents[name][0]].add_child(nodes[name])
        else:
            # identify this node as the root
            tree = nodes[name] 
            orphan_nodes += 1
    # there can only be one root
    assert orphan_nodes == 1
    # if we have a germline outgroup, reroot on it
    if args.naive is not None:
        # tree must have been outgroup rooted on the naive sequence
        assert nodes[args.naive] in tree.children
        tree.remove_child(nodes[args.naive])
        nodes[args.naive].add_child(tree)
        tree.dist = nodes[args.naive].dist
        nodes[args.naive].dist = 0
        tree = nodes[args.naive]

    if args.seed is not None:
        # need code here to color path from naive to seed and extract sequences along that sequence

    tree.link_to_alignment(alignment=aln, alg_format='fasta')
    tree.render(args.outfile+'.tree.png')

if __name__ == "__main__":
    main()

