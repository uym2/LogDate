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
from logdate.lca_lib import find_LCAs

MAX_ITER = 50000
MIN_NU = 1e-12
MIN_MU = 1e-5
EPSILON_t = 1e-5

def f_wLogDate(pseudo=0,seqLen=1000):
    def f(x,*args):
        return sum([sqrt(b+pseudo/seqLen)*log(abs(y))**2 for (y,b) in zip(x[:-2],args[0])])

    def g(x,*args):
        return np.array([2*sqrt(b+pseudo/seqLen)*log(abs(z))/z for (z,b) in zip(x[:-2],args[0])] + [0,0])

    def h(x,*args):
        return diags([sqrt(b+pseudo/seqLen)*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-2],args[0])]+[0,0])	

    return f,g,h


def setup_constraint(tree,smpl_times):
    active_set = [node for node in tree.postorder_node_iter() if node.is_active]
    N = len(active_set)-1
    cons_eq = []
    
    idx = 0
    b = [1.]*N
    
    for node in active_set:
        node.idx = idx
        idx += 1
        new_constraint = None        
        lb = node.taxon.label if node.is_leaf() else node.label
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
    
    return cons_eq,b

def logIt(tree,f_obj,cons_eq,b,x0=None,maxIter=MAX_ITER,pseudo=0,seqLen=1000,verbose=False):
    N = len([node for node in tree.postorder_node_iter() if node.is_active])-1

    bounds = Bounds(np.array([MIN_NU]*N+[MIN_MU]+[-np.inf]),np.array([np.inf]*(N+2)),keep_feasible=False)
    x_init = x0
    
    args = (b)
    linear_constraint = LinearConstraint(csr_matrix(cons_eq),[0]*len(cons_eq),[0]*len(cons_eq),keep_feasible=False)

    f,g,h = f_obj(pseudo=pseudo,seqLen=seqLen)
    
    print("Initial state:" )
    print("mu = " + str(x_init[-2]))
    print("fx = " + str(f(x_init,args)))
    #print(x_init)
    print("Maximum constraint violation: " + str(np.max(csr_matrix(cons_eq).dot(x_init))))  
    
    result = minimize(fun=f,method="trust-constr",x0=x_init,bounds=bounds,args=args,constraints=[linear_constraint],options={'disp': True,'verbose':3 if verbose else 1,'maxiter':maxIter},jac=g,hess=h)
   
    x_opt = result.x
    mu = x_opt[N]
    fx = result.fun  
    
    return mu,fx,x_opt  

def setup_smpl_time(tree,sampling_time=None,bw_time=False,root_time=0,leaf_time=1):
    smpl_times = {}   
    if root_time is not None: 
        smpl_times[tree.seed_node.label] = root_time if not bw_time else -root_time
    if leaf_time is not None:    
        for node in tree.leaf_nodes():
            smpl_times[node.taxon.label] = leaf_time if not bw_time else -leaf_time
    
    # case 1: no sampling time given --> return the smpl_times defined by root_time and leaf_time
    if not sampling_time:
        return smpl_times
    
    # case 2: read in user-specified sampling times; allow overriding root and leaf times
    queries = []
    times = []
    names = []
    with open(sampling_time,"r") as fin:
        for line in fin:
            ID,t = line.split()
            spl = ID.split('=')
            if len(spl) > 1:
                name,q = spl
            else:
                name = []
                q = spl[0]    
            q = q.split('+')
            t = float(t) if not bw_time else -float(t)
            queries.append(q)
            times.append(t)
            names.append(name)
    calibs = find_LCAs(tree,queries) 
    nodeIdx = 0
    for node,time,name in zip(calibs,times,names):
        if name:
            if node.is_leaf():
                node.taxon.label = name
            else:
                node.label = name
            lb = name
        else:
            lb = node.taxon.label if node.is_leaf() else node.label            
        smpl_times[lb] = time        
    return smpl_times   
    
def random_timetree(tree,sampling_time,nrep,seed=None,root_age=None,leaf_age=None,min_nleaf=3,fout=stdout):
    smpl_times = setup_smpl_time(tree,sampling_time=sampling_time,root_age=root_age,leaf_age=leaf_age)
    
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
    

def logDate_with_random_init(tree,f_obj,sampling_time=None,bw_time=False,root_time=0,leaf_time=1,nrep=1,min_nleaf=3,maxIter=MAX_ITER,seed=None,pseudo=0,seqLen=1000,verbose=False):
    smpl_times = setup_smpl_time(tree,sampling_time=sampling_time,bw_time=bw_time,root_time=root_time,leaf_time=leaf_time)    
    X,seed,T0 = random_date_init(tree,smpl_times,nrep,min_nleaf=min_nleaf,seed=seed)
    
    print("Finished initialization with random seed " + str(seed))
    f_min = None
    x_best = None

    i = 0
    n_succeed = 0
    
    cons_eq,b = setup_constraint(tree,smpl_times)

    for i,y in enumerate(zip(X,T0)):
        x0 = y[0] + [y[1]]
        _,f,x = logIt(tree,f_obj,cons_eq,b,x0=x0,maxIter=maxIter,pseudo=pseudo,seqLen=seqLen,verbose=verbose)
        print("Found local optimal for Initial point " + str(i+1))
        n_succeed += 1                
        
        if f_min is None or f < f_min:
            f_min = f
            x_best = x
            s_tree,t_tree = scale_tree(tree,x_best)
            compute_divergence_time(t_tree,smpl_times,bw_time=bw_time)
            print("Found a better log-scored configuration")
            print("New mutation rate: " + str(x_best[-2]))
            print("New log score: " + str(f_min))
    
    mu = x_best[-2]
    return mu,f_min,x_best,s_tree,t_tree 
    

def logDate_with_lsd(tree,sampling_time,root_age=None,brScale=False,lsdDir=None,seqLen=1000,maxIter=MAX_ITER):
    wdir = run_lsd(tree,sampling_time,outputDir=lsdDir)
    
    x0 = read_lsd_results(wdir)
    x1 = [1.]*len(x0)
    x1[-1] = x0[-1] 
    #x1=x0
    smpl_times = {}

    with open(sampling_time,"r") as fin:
        fin.readline()
        for line in fin:
            name,time = line.split()
            smpl_times[name] = float(time)


    #mu,f,x = calibrate_log_opt(tree,smpl_times,root_age=root_age,brScale=brScale,x0=x1)
    mu,f,x = logIt(tree,smpl_times,root_age=root_age,brScale=brScale,x0=x1,seqLen=seqLen,maxIter=maxIter)
    
    if lsdDir is None:
        rmtree(wdir)

    s_tree,t_tree = scale_tree(tree,x) 

    return mu,f,x,s_tree,t_tree   

def run_lsd(tree,sampling_time,outputDir=None):
    wdir = outputDir if outputDir is not None else mkdtemp()
    treefile = normpath(join(wdir,"mytree.tre"))
    tree_as_newick(tree,outfile=treefile,append=False)
    call([lsd_exec,"-i",treefile,"-d",sampling_time,"-v","-c"])
    return wdir
        

def read_lsd_results(inputDir):
# suppose LSD was run on the "mytree.newick" and all the outputs are placed inside inputDir
    log_file = normpath(join(inputDir, "mytree.tre.result")) 
    input_tree_file = normpath(join(inputDir, "mytree.tre")) 
    result_tree_file = normpath(join(inputDir, "mytree.tre.result.newick")) 

    s = open(log_file,'r').read()
    i = s.find("Tree 1 rate ") + 12
    mu = ""
    found_dot = False

    while (s[i] == '.' and not found_dot) or  (s[i] in [str(x) for x in range(10)]):
        mu += s[i]
        if s[i] == '.':
            found_dot = True
        i += 1
    mu = float(mu)

    taxa = TaxonNamespace()
    tree = Tree.get_from_path(input_tree_file,schema="newick",taxon_namespace=taxa,rooting="force-rooted") 
    tree.encode_bipartitions()
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2
    x0 = [10**-10]*N + [mu]
    
    idx = 0
    brlen_map = {}
    
    for node in tree.postorder_node_iter():
        if not node is tree.seed_node:
            key = node.bipartition
            brlen_map[key] = (idx,node.edge_length)
            idx += 1

    tree2 = Tree.get_from_path(result_tree_file,schema="newick",taxon_namespace=taxa,rooting="force-rooted")
    tree2.encode_bipartitions()
    
    for node in tree2.postorder_node_iter():
        if not node is tree2.seed_node:
            key = node.bipartition
            idx,el = brlen_map[key]
            if el > 0 and node.edge_length>0:
                x0[idx] = node.edge_length/float(el)

    return x0        
        

def calibrate_log_opt(tree,smpl_times,root_age=None,brScale=False,x0=None):
    def f0(x,*args):
        return sum([b*(w-1)*(w-1) for (w,b) in zip(x[:-1],args[0])])

    def f1(x):
        return sum([log(y)*log(y) for y in x[:-1]])
    
    def f2(x,*args):
        return sum([0.5*log(1+sqrt(b))*log(y)*log(y) for (y,b) in zip(x[:-1],args[0])])

    def g(x,a):    
        return sum([ x[i]*a[i] for i in range(len(x)) ])

    def h(x,p):
        return x[p]

    def g_gradient(x,a):
        return a

    def f1_gradient(x):
        return np.array([2*log(z)/z for z in x[:-1]] + [0])
    
    def f2_gradient(x,*args):
        return np.array([0.5*2*log(1+sqrt(b))*log(z)/z for (z,b) in zip(x[:-1],args[0])] + [0])
    
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2
    cons_eq = []
    cons_ineq = []
    
    idx = 0
    b = [1.]*N
    
    for i in range(N):
        cons_ineq.append({'type':'ineq','fun':h,'args':(i,)})

    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1
        if node.is_leaf():
            node.constraint = [0.0]*(N+1)
            node.constraint[node.idx] = node.edge_length
            node.constraint[N] = -smpl_times[node.taxon.label]
            b[node.idx] = node.edge_length
        else:
            children = list(node.child_node_iter())
                       
            a = np.array([ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N+1) ])
            cons_eq.append({'type':'eq','fun':g,'args':(a,)})

            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length
                b[node.idx] = node.edge_length
            elif root_age is not None:
                a = np.array(children[0].constraint[:-1] + [children[0].constraint[-1]-root_age])
                cons_eq.append({'type':'eq','fun':g,'args':(a,),'jac':g_gradient})    

    x0 = ([1.]*N + [0.01]) if x0 is None else x0
    bounds = [(0.00000001,999999)]*(N+1)
    args = (b)
    
    if brScale:
        result = minimize(fun=f2,x0=x0,args=args,bounds=bounds,constraints=cons_eq,method="SLSQP",options={'maxiter': 1500, 'ftol': 1e-10,'disp': True,'iprint':3},jac=f2_gradient)
    else:
        result = minimize(fun=f1,x0=x0,bounds=bounds,constraints=cons_eq,method="SLSQP",jac=f1_gradient ,options={'maxiter': 1500, 'ftol': 1e-10,'disp': True,'iprint':3})
    x = result.x
    s = x[N]
    
    
    #for node in tree.postorder_node_iter():
    #    if node is not tree.seed_node:
    #        node.edge_length *= x[node.idx]/s

    return s,f1(x),x

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
    
def compute_divergence_time(tree,sampling_time,bw_time=False):
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
        lb += "[t=" + (str(node.time) if not bw_time else str(-node.time)) + "]" 
        if not node.is_leaf():
            node.label = lb            
