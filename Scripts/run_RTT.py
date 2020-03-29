#! /usr/bin/env python

# implementation of the root-to-tip (RTT) method to compute t0 and mu

from sys import argv
import numpy as np
from dendropy import Tree

def compute_RTT_b(tree):
    b_map = {}
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            node.b = 0
        else:    
            node.b = node.parent_node.b + node.edge_length
        if node.is_leaf():
            b_map[node.taxon.label] = node.b
    return b_map
    
def run_RTT(tree,sampling_time):
    b_map = compute_RTT_b(tree)
    b = []
    t = []
    for x in sampling_time:
        b.append(b_map[x])
        t.append(sampling_time[x])

    b0,b1 = estimate_coef(np.array(t),np.array(b))
    mu = b1
    t0 = -b0/mu

    return t0,mu

def estimate_coef(x, y):
# do linear regression; copied from https://www.geeksforgeeks.org/linear-regression-python-implementation/
    # number of observations/points
    n = np.size(x)

    # mean of x and y vector
    m_x, m_y = np.mean(x), np.mean(y)

    # calculating cross-deviation and deviation about x
    SS_xy = np.sum(y*x) - n*m_y*m_x
    SS_xx = np.sum(x*x) - n*m_x*m_x

    # calculating regression coefficients
    b_1 = SS_xy / SS_xx
    b_0 = m_y - b_1*m_x

    return(b_0, b_1)

treefile = argv[1]  
timefile = argv[2]

smpl_time = {}
with open(timefile,'r') as fin:
    for line in fin:
        name,time = line.strip().split()
        smpl_time[name] = float(time)

with open(treefile,'r') as fin:
    for line in fin:
        tree = Tree.get(data=line,schema="newick",preserve_underscores=True)
        t0,mu = run_RTT(tree,smpl_time)
        print(t0,mu)
