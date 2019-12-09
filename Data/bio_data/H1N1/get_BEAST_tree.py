#! /usr/bin/env python

import treeswift
from sys import argv

infile = argv[1]
outfile = argv[2]


tree = treeswift.read_tree_nexus(infile)['TREE1']

for node in tree.traverse_preorder():
    lb = node.get_label().split('[')[0]
    node.set_label(lb)

tree.write_tree_newick(outfile)
