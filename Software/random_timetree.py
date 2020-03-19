#! /usr/bin/env python

import logdate
from logdate.logD_lib import random_timetree
from dendropy import Tree
import dendropy
#import treeswift
from logdate.tree_lib import tree_as_newick
import argparse
from sys import argv

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input tree")
parser.add_argument("-t","--samplingTime",required=False,help="Sampling time at leaf nodes. Default: None")
parser.add_argument("-p","--rep",required=False,help="The number of random replicates. Default: 1")
args = vars(parser.parse_args())

tree = Tree.get_from_path(args["input"],'newick',preserve_underscores=True)
sampling_time = args["samplingTime"]
nrep = int(args["rep"]) if args["rep"] else 1
#randseed = int(args["rseed"]) if args["rseed"] else None

random_timetree(tree,sampling_time,nrep)
