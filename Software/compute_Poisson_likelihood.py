#! /usr/bin/env


from dendropy import Tree, TaxonNamespace
from sys import argv
import argparse
from scipy.stats import poisson
from math import log

parser = argparse.ArgumentParser()

parser.add_argument("-b","--bTree", required=True,help="Input tree with branch length in substitution unit")
parser.add_argument("-c","--cTree", required=True,help="Calibrated tree. Must have the same rooted topology as bTree")
parser.add_argument("-r","--rate", required=False,help="Global mutation rate. Default: 1.0")
parser.add_argument("-l","--length", required=False,help="Sequence length. Default: 1000")
parser.add_argument("-o","--output", required=True,help="Output file")

args = vars(parser.parse_args())

taxa = TaxonNamespace()

b_tree = Tree.get_from_path(args["bTree"],schema="newick",taxon_namespace=taxa)
c_tree = Tree.get_from_path(args["cTree"],schema="newick",taxon_namespace=taxa)
r = float(args["rate"]) if args["rate"] else 1.0
l = float(args["length"]) if args["length"] else 1000

b_tree.is_rooted = True
c_tree.is_rooted = True

b_tree.encode_bipartitions()
c_tree.encode_bipartitions()

mapping = {}
llh = 0

for node in b_tree.postorder_node_iter():
    if node is not b_tree.seed_node:
        mapping[node.bipartition] = node.edge_length
        
        
for node in c_tree.postorder_node_iter():
    if node is not c_tree.seed_node:
        b = int(mapping[node.bipartition]*l+1)
        ld = int(node.edge_length*r*l+1)
        llh += log(poisson.pmf(b,ld))

with open(args["output"],'w') as fout:
    fout.write("Log-likelihood: " + str(llh))


