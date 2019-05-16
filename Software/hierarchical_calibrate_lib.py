#! /usr/bin/env python

from dendropy import Tree
import numpy as np
#from QP import quadprog_solve_qp, cvxopt_solve_qp
from math import exp,log, sqrt
from scipy.optimize import minimize

def calibrate_node(node,sampling_times):
    def f1(x,*args):
        w1,w2,mu = x
        return log(w1)*log(w1) + log(w2)*log(w2) + args[0]*log(mu)*log(mu) + args[1]*log(mu)
    
    def g1(x,*args):
        w1,w2,mu = x
        return args[0]*w1 - args[1]*w2 + args[2]*mu    

    def f2(x,*args):
        w1,w2,alpha1,alpha2 = x
        return log(w1)*log(w1) + log(w2)*log(w2) + args[0]*log(alpha1)*log(alpha1) + args[1]*log(alpha2)*log(alpha2) + 2*args[2]*log(alpha1) + 2*args[3]*log(alpha2)
    
    def g2(x,*args):
        w1,w2,alpha1,alpha2 = x
        return args[0]*w1 - args[1]*w2 + args[2]*alpha1 - args[3]*alpha2
    
    def f3(x,*args):
        w1,w2,alpha1,alpha2,mu = x
        return log(w1)*log(w1) + log(w2)*log(w2) + args[0]*log(alpha1)*log(alpha1) + args[1]*log(alpha2)*log(alpha2) + 2*args[2]*log(alpha1) + 2*args[3]*log(alpha2)
    
    def g3(x,*args):
        w1,w2,alpha1,alpha2,mu = x
        return args[0]*w1 - args[1]*w2 + args[2]*alpha1 - args[3]*alpha2 - args[4]*mu

    def g0(x):
        mu = x[-1]  
        return mu  

    if node.is_leaf():
        node.N = 0
        node.LG = 0
        node.h = 0
        node.t = sampling_times[node.taxon.label]
        node.mu = 1.0
        node.w = 1.0
    else:
        node1,node2 = node.child_nodes() # assuming node has exactly two children (the tree is perfectly bifurcating)
        if node1.mu is not None and node2.mu is not None:
            P = node1.h/node1.mu - node2.h/node2.mu - node1.t + node2.t
            Q = node1.LG + node2.LG - node1.N*node1.mu - node2.N*node2.mu
            args = (node1.N + node2.N, 2*Q)
            x0 = [1.0,1.0,node1.mu]
            bounds = [(0.00000001,999999)]*3
            w1,w2,mu = minimize(args=args,fun=f1,x0=x0,bounds=bounds,constraints=[{'type':'eq','fun':g1,'args':(node1.edge_length,node2.edge_length,P,)}],method="SLSQP").x 

            print("Scaling: " + str(node1.mu) + " " + str(node2.mu) + " " + str(mu))
        
            alpha1 = mu/node1.mu
            alpha2 = mu/node2.mu

        else:
            if node1.mu is None:
                # swap node1 and node2
                temp = node1
                node1 = node2
                node2 = temp
             # now we know node2.mu must be None
            if node1.mu is not None:
                args = (node1.N,node2.N,node1.LG,node2.LG)
                x0 = [1.0,1.0,1.0,1.0]
                bounds = [(0.00000001,999999)]*4
                
                w1,w2,alpha1,alpha2 = minimize(args=args,fun=f2,x0=x0,bounds=bounds,constraints=[{'type':'eq','fun':g2,'args':(node1.edge_length,node2.edge_length,node1.h-(node1.t-node2.t)*node1.mu,node2.h,)}],method="SLSQP").x                 
                mu = node1.mu*alpha1
            else:
                # both are None
                if node1.t == node2.t:
                    args = (node1.N,node2.N,node1.LG,node2.LG)
                    x0 = [1.0,1.0,1.0,1.0]
                    bounds = [(0.00000001,999999)]*4
                    w1,w2,alpha1,alpha2 = minimize(args=args,fun=f2,x0=x0,bounds=bounds,constraints=[{'type':'eq','fun':g2,'args':(node1.edge_length,node2.edge_length,node1.h,node2.h,)}],method="SLSQP").x 
                    mu = 1.0
                else:
                    args = (node1.N,node2.N,node1.LG,node2.LG)
                    x0 = [1.0,1.0,1.0,1.0,0.006]
                    bounds = [(0.00000001,999999)]*5
                    w1,w2,alpha1,alpha2,mu = minimize(args=args,fun=f3,x0=x0,bounds=bounds,constraints=[{'type':'eq','fun':g3,'args':(node1.edge_length,node2.edge_length,node1.h,node2.h,node1.t-node2.t,)}],method="SLSQP").x 

        node1.alpha = alpha1
        node2.alpha = alpha2
        node.N = node1.N + node2.N + 2
        node.LG = node1.N*log(node1.alpha) + node2.N*log(node2.alpha) + node1.LG + node2.LG + log(w1) + log(w2)
        node.h = node1.h*node1.alpha + w1*node1.edge_length
        node.t = node1.t
        node.mu = mu
        node1.w = w1
        node2.w = w2

def calibrate_bUp(a_tree,sampling_times):
    for node in a_tree.postorder_node_iter():
        calibrate_node(node,sampling_times)

def calibrate_tDown(a_tree):
    a_tree.seed_node.alpha = 1
    mu = a_tree.seed_node.mu
    print(mu)

    for node in a_tree.preorder_node_iter():
        if node is not a_tree.seed_node:
            node.alpha *= node.parent_node.alpha
            node.edge_length *= node.w*node.parent_node.alpha
        node.h *= node.alpha/mu

def compute_age(a_tree):
    for node in a_tree.postorder_node_iter():
        if node.is_leaf():
            node.age = node.t
        else:
            node.age = node.t - node.h 
            
def hierachical_logDate(a_tree,sampling_times): 
    calibrate_bUp(a_tree,sampling_times)
    calibrate_tDown(a_tree)
    compute_age(a_tree)
    
from sys import argv
from dendropy import Tree

a_tree = Tree.get_from_path(argv[1],"newick")
sampling_times = {}

with open(argv[2],'read') as fin:
    fin.readline()
    for line in fin:
        taxon,time = line.rstrip().split()
        sampling_times[taxon] = float(time)

hierachical_logDate(a_tree,sampling_times)              

#ID = 0
for node in a_tree.preorder_node_iter():
    #node.label = ID
    #ID += 1
    if node.is_leaf():
        print(str(node.taxon.label) + " " + str(node.age))
    else:
        print(str(node.label) + " " + str(node.age))
            
    
a_tree.write_to_path("Test.tre","newick")        
