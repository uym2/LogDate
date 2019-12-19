#! /usr/bin/env python

from dendropy import Tree
from sys import argv,stdout
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input time tree")
parser.add_argument("-m","--mode",required=True,help="Mode: either all_internal or root_age. If all_interal is chosen, the input tree MUST have distinct labels for ALL internal nodes")
parser.add_argument("-t","--samplingTime",required=True,help="Sampling time")
parser.add_argument("-o","--output",required=False,help="Output file")

args = vars(parser.parse_args())

tree = Tree.get_from_path(args["input"],'newick',preserve_underscores=True)
sampling_time = args["samplingTime"]
n = 0
smplt = {}
EPSILON = 1e-3
ages = {}
root_age = None

with open(sampling_time,'r') as fin:
    n = int(fin.readline())
    for line in fin:
        taxon,time = line.strip().split()
        assert taxon not in smplt, "repeat sampling time for taxon " + taxon
        smplt[taxon] = float(time)

#assert len(list(tree.leaf_node_iter())) == n, "number of taxa and sampling times do not match!"

for node in tree.postorder_node_iter():
    if node.is_leaf():
        assert node.taxon.label in smplt, "missing sampling time for taxon " + node.taxon.label
        node.age = smplt[node.taxon.label]
    else:
        age = None
        for ch in node.child_node_iter():
            if age is not None:
                assert abs(age - (ch.age - ch.edge_length)) < EPSILON, "inconsistent between sampling time and the time tree for taxon " + node.label
            else:
                age = ch.age - ch.edge_length
        node.age = age        
        if args["mode"] == "all_internal":
            ages[node.label] = node.age
        
        if node is tree.seed_node:
            root_age = node.age

fout = open(args["output"],'w') if args["output"] else stdout

if args["mode"] == "all_internal":
    for label in sorted(ages):
        fout.write(label + "\t" + str(ages[label]) + "\n")
else:
    fout.write(str(root_age))

if args["output"]:
    fout.close()   
