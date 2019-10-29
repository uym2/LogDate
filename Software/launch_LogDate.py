#! /usr/bin/env python

from logdate.logD_lib import calibrate_log_opt, read_lsd_results, logDate_with_lsd, logDate_with_random_init
from logdate.logD_CI_lib import logCI_with_lsd
from dendropy import TreeList
import treeswift
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-t","--samplingTime",required=False,help="Sampling time at leaf nodes. Default: None")
parser.add_argument("-r","--rootAge",required=False,help="Root age. Can be used with -t but by default -t will set root age to 0. Default: None if -t is specified else 0")
parser.add_argument("-f","--leafAge",required=False,help="Leaf age. To be used with root age to infer relative time. Will be overried by -t if it is specified. Default: None if -t is specified else 1.")
parser.add_argument("-o","--output",required=True,help="The output trees with branch lengths in time unit")
parser.add_argument("-d","--tempdir",required=False,help="The output from lsd will be kept in the specified directory")
parser.add_argument("-b","--brScale",required=False,help="Per-branch weighting strategy. Options include: 'sqrt', 'log', 'lsd'. Default: No weighting")
parser.add_argument("-p","--rep",required=False, help="The number of random replicates for initialization. Default: use 1 initial point")
parser.add_argument("-s","--rseed",required=False, help="Random seed to generate starting tree initial points")
parser.add_argument("-l","--seqLen",required=False, help="The length of the sequences. Default: 1000")
parser.add_argument("-m","--maxIter",required=False, help="The maximum number of iterations for optimization. Default: 50000")

args = vars(parser.parse_args())

myTrees = TreeList.get_from_path(args["input"],'newick',preserve_underscores=True)
sampling_time = args["samplingTime"]
rootAge = float(args["rootAge"]) if args["rootAge"] else None
leafAge = float(args["leafAge"]) if args["leafAge"] else None
lsdDir = args["tempdir"] if args["tempdir"] else None
nrep = int(args["rep"]) if args["rep"] else 1
seqLen = int(args["seqLen"]) if args["seqLen"] else 1000
brScale = args["brScale"]
maxIter = int(args["maxIter"]) if args["maxIter"] else 50000
randseed = int(args["rseed"]) if args["rseed"] else None

with open(args["output"],"w") as fout:
    for i,tree in enumerate(myTrees):
        print("Dating tree " + str(i+1))
        #mu,f,x,s_tree,t_tree= logDate_with_random_init(tree,sampling_time=sampling_time,root_age=rootAge,leaf_age=leafAge,brScale=brScale,seqLen=seqLen,nrep=nrep,min_nleaf=10,maxIter=maxIter,seed=randseed)
        mu,f,x,s_tree,t_tree = logDate_with_lsd(tree,sampling_time,root_age=rootAge,brScale=brScale,lsdDir=None,seqLen=seqLen,maxIter=maxIter)
        t_tree_swift = treeswift.read_tree_dendropy(t_tree)
        fout.write(t_tree_swift.newick())
        print("Clock rate: " + str(mu))
        print("Log score: " + str(f))

#if args["scaledTree"]:
#s_tree.write_to_path(args["scaledTree"],"newick")

