#! /usr/bin/env python

import logdate
from logdate.logD_lib import calibrate_log_opt, read_lsd_results, logDate_with_lsd, logDate_with_random_init,logDate_with_penalize_llh
from logdate.logD_lib import f_LF,f_lsd,f_PL,run_LF_cvxpy,f_wlogDate_linear_scale
from logdate.logD_lib import f_wlogDate_sqrt_scale, f_logDate_sqrt_scale, f_logDate_sqrt_b, f_logDate
from dendropy import Tree
import dendropy
#import treeswift
from logdate.tree_lib import tree_as_newick
import argparse
from sys import argv

print("Launching " + logdate.PROGRAM_NAME + " version " + logdate.PROGRAM_VERSION)
print(logdate.PROGRAM_NAME + " was called as follow")
print(" ".join(argv))

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input tree")
parser.add_argument("-j","--obj",required=False,help="Objective function. Either LogDate or wLogDate. Default: wLogDate")
parser.add_argument("-t","--samplingTime",required=False,help="Sampling time at leaf nodes. Default: None")
parser.add_argument("-r","--rootAge",required=False,help="Root age. Can be used with either -f or -t, but not both. Default: None if -t is specified else 0")
parser.add_argument("-f","--leafAge",required=False,help="Leaf age. To be used with root age to infer relative time. Will be overried by -t if -t is specified. Default: None if -t is specified else 1.")
parser.add_argument("-o","--output",required=True,help="The output trees with branch lengths in time unit.")
parser.add_argument("-v","--verbose",action='store_true',help="Show verbose message. Default: NO")
parser.add_argument("-p","--rep",required=False,help="The number of random replicates for initialization. Default: use 1 initial point")
parser.add_argument("-s","--rseed",required=False,help="Random seed to generate starting tree initial points")
parser.add_argument("-l","--seqLen",required=False,help="The length of the sequences. Default: 1000")
parser.add_argument("-m","--maxIter",required=False,help="The maximum number of iterations for optimization. Default: 50000")
parser.add_argument("-u","--addpseudo",required=False,help="Add pseudo counting for per-branch weighting.Default: 0.01")
parser.add_argument("-z","--zero",required=False,help="Set zero-length branches (if any) to this number. LogDate cannot process zero-length branches. Default: 1e-10")

args = vars(parser.parse_args())

tree = Tree.get_from_path(args["input"],'newick',preserve_underscores=True)
sampling_time = args["samplingTime"]
rootAge = float(args["rootAge"]) if args["rootAge"] else None
leafAge = float(args["leafAge"]) if args["leafAge"] else None

if sampling_time is None and rootAge is None and leafAge is None:
    rootAge = 0
    leafAge = 1

nrep = int(args["rep"]) if args["rep"] else 1
seqLen = int(args["seqLen"]) if args["seqLen"] else 1000
pseudo = 0.01 if args["addpseudo"] is None else float(args["addpseudo"])
maxIter = int(args["maxIter"]) if args["maxIter"] else 50000
randseed = int(args["rseed"]) if args["rseed"] else None
zero_len = float(args["zero"]) if args["zero"] else 1e-10

brScale = None
f_obj = f_logDate if args["obj"] == "LogDate" else f_logDate_sqrt_b
objective = "wLogDate" if args["obj"] is None else args["obj"]
verbose = args["verbose"]

for node in tree.postorder_node_iter():
    if node is not tree.seed_node and node.edge_length == 0:
        node.edge_length = zero_len

mu,f,x,s_tree,t_tree = logDate_with_random_init(tree,f_obj,sampling_time=sampling_time,root_age=rootAge,leaf_age=leafAge,nrep=nrep,min_nleaf=10,maxIter=maxIter,seed=randseed,pseudo=pseudo,seqLen=seqLen,verbose=verbose)
tree_as_newick(t_tree,outfile=args["output"],append=False)

print("Clock rate: " + str(mu))
print("Log score: " + str(f))
