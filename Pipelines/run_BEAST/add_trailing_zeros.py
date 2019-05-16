#! /usr/bin/env python

from dendropy import Tree
from sys import argv

infile = argv[1]
outfile = argv[2]
ndigits = int(argv[3])

with open(infile,'r') as f:
    with open(outfile,'w') as fout:
        for line in f:
            t = Tree.get(data=line,schema="newick")
            for node in t.leaf_nodes():
                node.taxon.label = node.taxon.label.zfill(ndigits)
            fout.write(t.as_string("newick"))    
