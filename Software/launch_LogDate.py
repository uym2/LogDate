#! /usr/bin/env python

from logdate.logD_lib import calibrate_log_opt, read_lsd_results, logDate_with_lsd, logDate_with_random_init
from logdate.logD_CI_lib import logCI_with_lsd
from dendropy import TreeList
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-t","--samplingTime",required=True,help="Sampling time")
parser.add_argument("-r","--rootAge",required=False,help="Root age")
parser.add_argument("-o","--output",required=True,help="The output trees with branch lengths in time unit")
parser.add_argument("-s","--scaledTree",required=False,help="The output trees with branch lengths scaled")
parser.add_argument("-d","--tempdir",required=False,help="The output from lsd will be kept in the specified directory")
parser.add_argument("-c","--CI",required=False,action='store_true',help="Use confidence interval of Poisson in the objective function. Default: NO")
parser.add_argument("-b","--brScale",required=False,help="Per-branch weighting strategy. Options include: 'sqrt', 'log', 'lsd'. Default: No weighting")
parser.add_argument("-p","--rep",required=False, help="The number of random replicates for initialization. Default: use lsd initialization instead")
parser.add_argument("-l","--seqLen",required=False, help="The length of the sequences. Default: 1000")

args = vars(parser.parse_args())

myTrees = TreeList.get_from_path(args["input"],'newick')
sampling_time = args["samplingTime"]
rootAge = float(args["rootAge"]) if args["rootAge"] else None
lsdDir = args["tempdir"] if args["tempdir"] else None
nrep = int(args["rep"]) if args["rep"] else None
useCI = args["CI"]
seqLen = int(args["seqLen"]) if args["seqLen"] else 1000
brScale = args["brScale"]

for tree in myTrees:
    if useCI:
        mu,f,x,s_tree,t_tree = logCI_with_lsd(tree,sampling_time,root_age=rootAge,seqLen=seqLen,lsdDir=lsdDir)
    else:    
        if nrep is None:
            mu,f,x,s_tree,t_tree = logDate_with_lsd(tree,sampling_time,root_age=rootAge,brScale=brScale,lsdDir=lsdDir,seqLen=seqLen)
        else:
            mu,f,x,s_tree,t_tree = logDate_with_random_init(tree,sampling_time,root_age=rootAge,brScale=brScale,seqLen=seqLen,nrep=nrep,min_nleaf=10)

t_tree.write_to_path(args["output"],"newick")
if args["scaledTree"]:
    s_tree.write_to_path(args["scaledTree"],"newick")

print("Clock rate: " + str(mu))
print("Log score: " + str(f))
