#! /usr/bin/env


from treeswift import *
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input tree")
parser.add_argument("-a","--alpha",required=False,help="Scaling factor for each branch. Default: 1/h where h is the tree height")
parser.add_argument("-o","--output",required=True,help="Output tree")

args = vars(parser.parse_args())

infile = args["input"]
alpha = float(args["alpha"]) if args["alpha"] else None # multiply each branch by this number
outfile = args["output"]

tree = read_tree_newick(infile)

if alpha is None:
    # compute tree height h and set alpha to 1/h
    h = tree.height()
    alpha = 1/h

tree.scale_edges(alpha)
tree.write_tree_newick(outfile)
