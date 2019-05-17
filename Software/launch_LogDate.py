#! /usr/bin/env python

from logdate.logD_lib import calibrate_log_opt, read_lsd_results, logDate_with_lsd
from dendropy import TreeList
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-s","--samplingTime",required=True,help="Sampling time")
parser.add_argument("-r","--rootAge",required=False,help="Root age")
parser.add_argument("-o","--output",required=True,help="The output trees with branch lengths in time unit")
parser.add_argument("-d","--tempdir",required=False,help="The output from lsd will be kept in the specified directory")

args = vars(parser.parse_args())

myTrees = TreeList.get_from_path(args["input"],'newick')
sampling_time = args["samplingTime"]
rootAge = float(args["rootAge"]) if args["rootAge"] else None
lsdDir = args["tempdir"] if args["tempdir"] else None

for tree in myTrees:
    s,f = logDate_with_lsd(tree,sampling_time,root_age=rootAge,brScale=False,lsdDir=lsdDir)

myTrees.write_to_path(args["output"],"newick")
print("Clock rate: " + str(s))
print("Log score: " + str(f))
