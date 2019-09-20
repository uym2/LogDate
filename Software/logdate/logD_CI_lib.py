from dendropy import Tree,TaxonNamespace
import numpy as np
from math import exp,log, sqrt
from scipy.stats.stats import pearsonr
from scipy.stats import chi2
from scipy.optimize import minimize, LinearConstraint,Bounds
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
from subprocess import check_output,call
from tempfile import mkdtemp
from shutil import copyfile, rmtree
from os import remove
from copy import deepcopy
from logdate.init_lib import random_date_init
from logdate.logD_lib import run_lsd,read_lsd_results,scale_tree

def poisson_interval(k, alpha=0.05): 
    """
    uses chisquared info to get the poisson interval. Uses scipy.stats 
    (imports in function). 
    """
    a = alpha
    low, high = (chi2.ppf(a/2, 2*k) / 2, chi2.ppf(1-a/2, 2*k + 2) / 2)
    if k == 0: 
        low = 0.0
    return low, high

def setup_constr(tree,smpl_times,root_age=None,seqLen=1000):
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2
    cons_eq = []
    
    idx = 0
    B = [0]*N

    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1
        if node.is_leaf():
            node.constraint = [0.0]*(N+1)
            node.constraint[node.idx] = node.edge_length
            node.constraint[N] = -smpl_times[node.taxon.label]
            B[node.idx] = int(seqLen*node.edge_length)
        else:
            children = list(node.child_node_iter())           
            a = [ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N+1) ]
            cons_eq.append(a)

            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length
                B[node.idx] = int(seqLen*node.edge_length)
            elif root_age is not None:
                a = children[0].constraint[:-1] + [children[0].constraint[-1]-root_age]
                cons_eq.append(a)    
    
    bounds = Bounds(np.array([1e-3]*(N+1)),np.array([9999999]*(N+1)))
    linear_constraint = LinearConstraint(cons_eq,[0]*(n-1),[0]*(n-1))

    return bounds, linear_constraint, B

def logCI(bounds,linear_constraint,B,x0=None):
    def g(y,l,u):
        if y < l:
            return (log(abs(y))-log(abs(l)))**2;
        elif y > u:
            return (log(abs(y))-log(abs(u)))**2;
        else:
            return 0
            
    def g_first(y,l,u):
        if y < l:
            return 2*(log(y)-log(l))/y
        elif y > u:    
            return 2*(log(y)-log(u))/y
        else:
            return 0
    
    def f(x,*args):
        CI = args[0]
        return sum([ g(y,l,u) for (y,(l,u)) in zip(x[:-1],CI) ])
            
    def f_gradient(x,*args):
        CI = args[0]
        return np.array([ g_first(y,l,u) for (y,(l,u)) in zip(x[:-1],CI) + [0]  ])
            
    x0 = ([1.]*N + [0.01]) if x0 is None else x0
    CI = poisson_interval(B)
    args = (CI)

    result = minimize(fun=f,method="trust-constr",x0=x0,bounds=bounds,args=args,constraints=[linear_constraint],jac=f_gradient,options={'disp':True,'verbose':3,'maxiter':5000})
   
    x = result.x
    mu = x[-1]
   
    fx = f(x,args)
    return mu,fx,x                 

def logCI_with_lsd(tree,sampling_time,root_age=None,seqLen=1000,lsdDir=None):
    wdir = run_lsd(tree,sampling_time,outputDir=lsdDir)
    
    x0 = read_lsd_results(wdir)
    x1 = [1.]*len(x0)
    x1[-1] = x0[-1] 
    smpl_times = {}

    with open(sampling_time,"r") as fin:
        fin.readline()
        for line in fin:
            name,time = line.split()
            smpl_times[name] = float(time)
    

    bounds, linear_constraint, B = setup_constr(tree,smpl_times,root_age=root_age,seqLen=seqLen)
    mu,f,x = logCI(bounds,linear_constraint,B,x0=x1)
    
    if lsdDir is None:
        rmtree(wdir)

    s_tree,t_tree = scale_tree(tree,x) 

    return mu,f,x,s_tree,t_tree   
