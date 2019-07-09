#! /usr/bin/env python

from logdate.logD_lib import calibrate_log_opt, read_lsd_results, logDate_with_lsd, logDate_with_random_init
from dendropy import TreeList
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-t","--samplingTime",required=True,help="Sampling time")
parser.add_argument("-r","--rootAge",required=False,help="Root age")
parser.add_argument("-o","--output",required=True,help="The output trees with branch lengths in time unit")
parser.add_argument("-s","--scaledTree",required=False,help="The output trees with branch lengths scaled")
parser.add_argument("-d","--tempdir",required=False,help="The output from lsd will be kept in the specified directory")
parser.add_argument("-b","--brScale",required=False,action='store_true',help="Do branch length scaling. Default: NO")
parser.add_argument("-p","--rep",required=False, help="The number of random replicates for initialization. Default: use lsd initialization instead")

args = vars(parser.parse_args())

myTrees = TreeList.get_from_path(args["input"],'newick')
sampling_time = args["samplingTime"]
rootAge = float(args["rootAge"]) if args["rootAge"] else None
lsdDir = args["tempdir"] if args["tempdir"] else None
nrep = int(args["rep"]) if args["rep"] else None

for tree in myTrees:
    if nrep is None:
        mu,f,x,s_tree,t_tree = logDate_with_lsd(tree,sampling_time,root_age=rootAge,brScale=args["brScale"],lsdDir=lsdDir)
    else:
        mu,f,x,s_tree,t_tree = logDate_with_random_init(tree,sampling_time,root_age=rootAge,brScale=args["brScale"],nrep=rep,min_nleaf=10)

t_tree.write_to_path(args["output"],"newick")
if args["scaledTree"]:
    s_tree.write_to_path(args["scaledTree"],"newick")

print("Clock rate: " + str(mu))
print("Log score: " + str(f))
