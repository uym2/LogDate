#! /usr/bin/env python

from logdate.logD_lib import calibrate_log_opt
from dendropy import TreeList
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-s","--samplingTime",required=True,help="Sampling time")
parser.add_argument("-r","--rootAge",required=False,help="Root age")
parser.add_argument("-o","--output",required=True,help="The output trees with branch lengths in time unit")

args = vars(parser.parse_args())

myTrees = TreeList.get_from_path(args["input"],'newick')
smpl_times = {}
rootAge = float(args["rootAge"]) if args["rootAge"] else None

with open(args["samplingTime"],"r") as fin:
    fin.readline()
    for line in fin:
        name,time = line.split()
        smpl_times[name] = float(time)

for tree in myTrees:
    s = calibrate_log_opt(tree,smpl_times,root_age=rootAge,brScale=False)

myTrees.write_to_path(args["output"],"newick")
print("Clock rate: " + str(s))
