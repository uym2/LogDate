from dendropy import Tree,TaxonNamespace
import numpy as np
from math import exp,log, sqrt
from scipy.optimize import minimize, LinearConstraint,Bounds
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists
from subprocess import check_output,call
from tempfile import mkdtemp
from shutil import copyfile, rmtree
from os import remove
from copy import deepcopy
from logdate.fixed_init_lib import random_date_init
#from logdate.init_lib_old import random_date_init
import platform
from scipy.sparse import diags
from scipy.sparse import csr_matrix
#import cvxpy as cp
import dendropy
from logdate.tree_lib import tree_as_newick
from sys import stdout

MAX_ITER = 50000
MIN_NU = 1e-12
MIN_MU = 1e-5
EPSILON_t = 1e-5
EPSILON = 1e-5

def f_wLogDate(pseudo=0,seqLen=1000):
    def f(x,*args):
        return sum([sqrt(b+pseudo/seqLen)*log(abs(y))**2 for (y,b) in zip(x[:-2],args[0])])

    def g(x,*args):
        return np.array([2*sqrt(b+pseudo/seqLen)*log(abs(z))/z for (z,b) in zip(x[:-2],args[0])] + [0,0])

    def h(x,*args):
        return diags([sqrt(b+pseudo/seqLen)*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-2],args[0])]+[0,0])	

    return f,g,h


def add_ieq_constraint(tree,a_node,tmin,tmax,N):
    curr_node = a_node
    lower_constr = [0.0]*N + [-tmin,1] if tmin != -np.inf else None
    upper_constr = [0.0]*N + [-tmax,1] if tmax != np.inf else None

    while curr_node is not tree.seed_node:
        if lower_constr is not None:
            lower_constr[curr_node.idx] = curr_node.edge_length
        if upper_constr is not None:
            upper_constr[curr_node.idx] = curr_node.edge_length    
        curr_node = curr_node.parent_node

    return lower_constr,upper_constr

def setup_constraint(tree,smpl_times,soft_calibs):
# convention: smpl_times contains hard calibrations, soft_calibs contains soft calibrations    
    active_set = [node for node in tree.postorder_node_iter() if node.is_active]
    N = len(active_set)-1
    cons_eq = [] # equality constraints
    cons_lower = [] # inequality constraints
    cons_upper = [] # inequality constraints
        
    b = [1.]*N
    
    # index tree
    idx = 0
    for node in active_set:
        node.idx = idx
        idx += 1
    
    for node in active_set:
        lb = node.taxon.label if node.is_leaf() else node.label
        
        # add inequality constraints
        if lb in soft_calibs:
            tmin,tmax = soft_calibs[lb]
            lower_constr,upper_constr = add_ieq_constraint(tree,node,tmin,tmax,N)
            if lower_constr is not None:
                cons_lower.append(lower_constr)
            if upper_constr is not None:
                cons_upper.append(upper_constr)    

        # add equality constraints
        new_constraint = None        
        if lb in smpl_times:
            node.height = 0
            new_constraint = [0.0]*(N+2)
            new_constraint[N] = -smpl_times[lb]            
            new_constraint[N+1] = 1           
        
        if not node.as_leaf:
            C = [ c for c in node.child_node_iter() if c.is_active ]
            # assumming each node has either one child or two children
            if len(C) == 1: # if it has one child
                c1 = C[0]
                c2 = None
            else: # it has exactly two children
                c1,c2 = C

            f0 = lb in smpl_times
            f1 = c1.constraint is not None
            f2 = c2 is not None and c2.constraint is not None

            if f1 and f2:
                if not f0:
                    node.height = min(c1.height,c2.height) + 1
                    new_constraint = c1.constraint if c1.height < c2.height else c2.constraint
                    a = [ (c1.constraint[i] - c2.constraint[i]) for i in range(N+2) ]
                    cons_eq.append(a)
                else:    
                    a1 = c1.constraint[:-2] + [smpl_times[lb]+c1.constraint[-2]] + [0]
                    a2 = c2.constraint[:-2] + [smpl_times[lb]+c2.constraint[-2]] + [0]
                    cons_eq.append(a1)
                    cons_eq.append(a2)
            elif f1:
                if not f0:
                    node.height = c1.height
                    new_constraint = c1.constraint
                else:    
                    a1 = c1.constraint[:-2] + [smpl_times[lb]+c1.constraint[-2]] + [0]
                    cons_eq.append(a1)
            elif f2:
                if not f0:
                    node.height = c2.height
                    new_constraint = c2.constraint        
                else:        
                    a2 = c2.constraint[:-2] + [smpl_times[lb]+c2.constraint[-2]] + [0]
                    cons_eq.append(a2)   
            
        node.constraint = new_constraint
        if node is not tree.seed_node and node.constraint is not None:    
            node.constraint[node.idx] = node.edge_length
            b[node.idx] = node.edge_length
    return cons_eq,cons_lower,cons_upper,b

def logIt(tree,f_obj,cons_eq,cons_lower,cons_upper,b,x0=None,maxIter=MAX_ITER,pseudo=0,seqLen=1000,verbose=False):
    N = len([node for node in tree.postorder_node_iter() if node.is_active])-1

    bounds = Bounds(np.array([MIN_NU]*N+[MIN_MU]+[-np.inf]),np.array([np.inf]*(N+2)),keep_feasible=True)
    x_init = x0
    
    args = (b)
    neq = len(cons_eq)
    nlower = len(cons_lower)
    nupper = len(cons_upper)
    
    #l = [0]*(neq+nlower) + [-np.inf]*nupper
    #r = [0]*neq + [np.inf]*nlower +  [0]*nupper
    
    #l = [0]*neq
    #r = [0]*neq

    #linear_constraint = LinearConstraint(csr_matrix(cons_eq+cons_lower+cons_upper),l,r,keep_feasible=True)
    eq_constraint = LinearConstraint(csr_matrix(cons_eq),[0]*neq,[0]*neq,keep_feasible=False)
    lower_constraint = LinearConstraint(csr_matrix(cons_lower),[0]*nlower,[np.inf]*nlower,keep_feasible=False)
    upper_constraint = LinearConstraint(csr_matrix(cons_upper),[-np.inf]*nupper,[0]*nupper,keep_feasible=False)

    f,g,h = f_obj(pseudo=pseudo,seqLen=seqLen)
    
    print("Initial state:" )
    print("mu = " + str(x_init[-2]))
    print("fx = " + str(f(x_init,args)))
    #print(x_init)
    print("Maximum constraint violation: " + str(np.max(csr_matrix(cons_eq).dot(x_init))))  
    
    result = minimize(fun=f,method="trust-constr",x0=x_init,bounds=bounds,args=args,constraints=[eq_constraint,lower_constraint,upper_constraint],options={'disp': True,'verbose':3 if verbose else 1,'maxiter':maxIter},jac=g,hess=h)
   
    x_opt = result.x
    mu = x_opt[N]
    fx = result.fun  
    
    return mu,fx,x_opt  

def setup_smpl_time(tree,sampling_time):
    smpl_times = {}
    soft_calib = {}
    with open(sampling_time,"r") as fin:
        for line in fin:
            item = line.strip().split()
            if len(item) < 2:
                print("Warning: ignore one invalid sampling time")
            elif len(item) == 2:    
                name,time = item
                smpl_times[name] = float(time)
            elif len(item) == 3:
                name,tmin,tmax = item
                tmin = float(tmin) if tmin is not None else -np.inf
                tmax = float(tmax) if tmax is not None else np.inf
                soft_calib[name] = (tmin,tmax)  
                #smpl_times[name] = (tmin+tmax)/2 # TEMP 
    return smpl_times,soft_calib   
    
def random_timetree(tree,sampling_time,nrep,seed=None,root_age=None,leaf_age=None,min_nleaf=3,fout=stdout):
    smpl_times,soft_calib = setup_smpl_time(tree,sampling_time=sampling_time,root_age=root_age,leaf_age=leaf_age)
    
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            node.fixed_age = smpl_times[node.taxon.label]
        else:    
            node.fixed_age = None
    
    setup_constraint(tree,smpl_times,root_age=root_age)
    X,seed,_ = random_date_init(tree,smpl_times,nrep,min_nleaf=min_nleaf,rootAge=root_age,seed=seed)
    print("Finished initialization with random seed " + str(seed))
    
    for x in X:
        s_tree,t_tree = scale_tree(tree,x)
        fout.write(t_tree.as_string("newick"))
    

def logDate_with_random_init(tree,f_obj,sampling_time,nrep=1,min_nleaf=3,maxIter=MAX_ITER,seed=None,pseudo=0,seqLen=1000,verbose=False):
    smpl_times,soft_calib = setup_smpl_time(tree,sampling_time)
    
    
    fixed_points = {}
    for x in smpl_times:
        fixed_points[x] = smpl_times[x]
    
    # quick and dirty solution
    for x in soft_calib:
        tmin,tmax = soft_calib[x]
        t = (tmin+tmax)/2        
        fixed_points[x] = t        

    X,seed,T0 = random_date_init(tree,fixed_points,nrep,min_nleaf=min_nleaf,seed=seed)
    print("Finished initialization with random seed " + str(seed))
    
    f_min = None
    x_best = None
    cons_eq,cons_lower,cons_upper,b = setup_constraint(tree,smpl_times,soft_calib)

    for i,y in enumerate(zip(X,T0)):
        mu = y[0][-1]
        t0 = y[1]
        b0 = mu*t0
        x0 = y[0] + [b0]
        _,f,x = logIt(tree,f_obj,cons_eq,cons_lower,cons_upper,b,x0=x0,maxIter=maxIter,pseudo=pseudo,seqLen=seqLen,verbose=verbose)
        print("Found local optimal for Initial point " + str(i+1))
        
        if f_min is None or f < f_min:
            f_min = f
            x_best = x
            s_tree,t_tree = scale_tree(tree,x_best)
            compute_divergence_time(t_tree,smpl_times)
            print("Found a better log-scored configuration")
            print("New mutation rate: " + str(x_best[-2]))
            print("New log score: " + str(f_min))
            print("Time tree")
            print(t_tree.as_string("newick"))
    
    mu = x_best[-2]
    print(x_best[-1]/mu)
    return mu,f_min,x_best,s_tree,t_tree 

def scale_tree(tree,x):
    taxa = tree.taxon_namespace
    s_tree = Tree.get(data=tree.as_string("newick"),taxon_namespace=taxa,schema="newick",rooting="force-rooted")

    tree.is_rooted = True
    tree.encode_bipartitions()
    s_tree.encode_bipartitions()

    mapping = {}
    mu = x[-2]
    for node in tree.postorder_node_iter():
        if node is not tree.seed_node and node.is_active:
            key = node.bipartition    
            mapping[key] = node.idx

    for node in s_tree.postorder_node_iter():
        if node is not s_tree.seed_node:
            if node.bipartition in mapping:
                idx = mapping[node.bipartition]
                node.edge_length *= x[idx]

    #t_tree = Tree.get(data=s_tree.as_string("newick"),taxon_namespace=taxa,schema="newick",rooting="force-rooted")
    t_tree = Tree.get(data=s_tree.as_string("newick"),schema="newick",rooting="force-rooted")
    
    for node in t_tree.postorder_node_iter():
        if node is not t_tree.seed_node:
            node.edge_length /= mu

    return s_tree,t_tree    
    
def compute_divergence_time(tree,sampling_time):
# compute and place the divergence time onto the node label of the tree
# must have at least one sampling time. Assumming the tree branches have been
# converted to time unit and are consistent with the given sampling_time
    calibrated = []
    for node in tree.postorder_node_iter():
        node.time = None
        lb = node.taxon.label if node.is_leaf() else node.label
        if lb in sampling_time:
            node.time = sampling_time[lb]
            calibrated.append(node)

    stk = []
    # push to stk all the uncalibrated nodes that are linked to (i.e. is parent or child of) any node in the calibrated list
    for node in calibrated:
        p = node.parent_node
        if p is not None and p.time is None:
            stk.append(p)
        if not node.is_leaf():
            stk += [ c for c in node.child_node_iter() if c.time is None ]            
    
    # compute divergence time of the remaining nodes
    while stk:
        node = stk.pop()
        lb = node.taxon.label if node.is_leaf() else node.label
        p = node.parent_node
        t = None
        if p is not None:
            if p.time is not None:
                t = p.time + node.edge_length
            else:
                stk.append(p)    
        for c in node.child_node_iter():
            if c.time is not None:
                t1 = c.time - c.edge_length
                t = t1 if t is None else t
                assert abs(t-t1) < EPSILON_t, "Inconsistent divergence time computed for node " + lb
            else:
                stk.append(c)
        node.time = t            
        
    # place the divergence time onto the label
    for node in tree.postorder_node_iter():                
        lb = node.taxon.label if node.is_leaf() else node.label
        assert node.time is not None, "Failed to compute divergence time for node " + lb 
        lb += "=" + str(node.time)
        if node.is_leaf():
            node.taxon.label = lb
        else:
            node.label = lb    
                            
