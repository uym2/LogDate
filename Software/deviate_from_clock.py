#! /usr/bin/env python

from calibration_lib import gamma_deviation_from_clock, exp_deviation_from_clock,lnorm_deviation_from_clock
from sys import argv
from dendropy import Tree

intrees = argv[1]
shape = argv[2] # e for exponential, g[alpha] for gamma with the specified alpha
rate = float(argv[3]) 
outtrees = argv[4]

outs = []
with open(intrees,'r') as fin:
    for line in fin:
        tree = Tree.get(data=line,schema="newick")
        if shape[0] == 'e':
            tree1 = exp_deviation_from_clock(tree,rate)
        elif shape[0] == 'g':
            alpha = float(shape[1:])    
            tree1 = gamma_deviation_from_clock(tree,alpha,rate)
        elif shape[0] == 'l':
            sd = float(shape[1:])
            tree1 = lnorm_deviation_from_clock(tree,sd,rate)
                
        outs.append(tree1.as_string("newick"))

with open(outtrees,'w') as fout:
    for tree in outs:
        fout.write(tree)
