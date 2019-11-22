#! /usr/bin/env python

import logdate
from logdate.logD_lib import calibrate_log_opt, read_lsd_results, logDate_with_lsd, logDate_with_random_init, f_LF, f_lsd
from logdate.logD_CI_lib import logCI_with_lsd
from logdate.logD_extend_lib import write_time_tree,log_from_random_init
from dendropy import Tree
import treeswift
import argparse

print("Launching " + logdate.PROGRAM_NAME + " version " + logdate.PROGRAM_VERSION)

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input tree")
parser.add_argument("-j","--obj",required=False,help="Objective function. Either LogDate, wLogDate, LF, or LSD. Default: wLogDate")
parser.add_argument("-t","--samplingTime",required=False,help="Sampling time at leaf nodes. Default: None")
parser.add_argument("-r","--rootAge",required=False,help="Root age. Can be used with either -f or -t, but not both. Default: None if -t is specified else 0")
parser.add_argument("-f","--leafAge",required=False,help="Leaf age. To be used with root age to infer relative time. Will be overried by -t if -t is specified. Default: None if -t is specified else 1.")
parser.add_argument("-o","--output",required=True,help="The output trees with branch lengths in time unit.")
parser.add_argument("-d","--tempdir",required=False,help="The output from lsd will be kept in the specified directory")
parser.add_argument("-e","--pseudo",action='store_true',help="Can only be used with wLogDate. Control pseudo-count branches. Default: NO")
parser.add_argument("-p","--rep",required=False,help="The number of random replicates for initialization. Default: use 1 initial point")
parser.add_argument("-s","--rseed",required=False,help="Random seed to generate starting tree initial points")
parser.add_argument("-l","--seqLen",required=False,help="The length of the sequences. Default: 1000")
parser.add_argument("-m","--maxIter",required=False,help="The maximum number of iterations for optimization. Default: 50000")

args = vars(parser.parse_args())

tree = Tree.get_from_path(args["input"],'newick',preserve_underscores=True)
sampling_time = args["samplingTime"]
rootAge = float(args["rootAge"]) if args["rootAge"] else None
leafAge = float(args["leafAge"]) if args["leafAge"] else None
lsdDir = args["tempdir"] if args["tempdir"] else None
nrep = int(args["rep"]) if args["rep"] else 1
seqLen = int(args["seqLen"]) if args["seqLen"] else 1000
maxIter = int(args["maxIter"]) if args["maxIter"] else 50000
randseed = int(args["rseed"]) if args["rseed"] else None

brScale = None
f_obj = None

if args["obj"] == "LF":
    f_obj = f_LF
elif args["obj"] == "LSD":
    f_obj = f_lsd
elif args["obj"] != "LogDate":
    brScale = "sqrt"            

if args["pseudo"]:
    smpl_time = {}
    with open(sampling_time,'r') as fin:
        fin.readline()
        for line in fin:
            taxon,time = line.split()
            smpl_time[taxon] = float(time)
    x,f = log_from_random_init(tree,smpl_time,root_age=rootAge,leaf_age=leafAge,brScale=brScale,nrep=nrep,min_nleaf=10,maxIter=maxIter,seed=randseed)
    write_time_tree(tree,outfile=args["output"])
    print("Best log-scored solution: ")
    print("Log score: " + str(f))
    print("Clock rate: " + str(x[1]))
    print("Root age: " + str(x[0]/x[1]))
else:    
    mu,f,x,s_tree,t_tree= logDate_with_random_init(tree,sampling_time=sampling_time,root_age=rootAge,leaf_age=leafAge,brScale=brScale,seqLen=seqLen,nrep=nrep,min_nleaf=3,maxIter=maxIter,seed=randseed,f_obj=f_obj)
    t_tree_swift = treeswift.read_tree_dendropy(t_tree)
    t_tree_swift.write_tree_newick(args["output"])
    print("Clock rate: " + str(mu))
    print("Log score: " + str(f))
