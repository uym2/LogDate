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
from logdate.init_lib import random_date_init
import platform
from scipy.sparse import diags
from scipy.sparse import csr_matrix
from logdate.util_lib import *
import treeswift

MAX_ITER = 50000
MIN_RATE = 1e-5
MIN_NU = 1e-10

lsd_file = "../lsd-0.2/bin/lsd.exe" if platform.system() == "Linux" else "../lsd-0.2/src/lsd"

lsd_exec=normpath(join(dirname(realpath(__file__)),lsd_file))

def f_logDate_lsd(c=10,s=1000):
# follow the LSD paper for weighting strategy
    def f(x,*args):
        return sum([(b+c/s)/s*log(abs(y))**2 for (y,b) in zip(x[:-1],args[0])])
    
    def g(x,*args):
        return np.array([2*(b+c/s)/s*log(abs(z))/z for (z,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        #return np.diag([(b+c/s)/s*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])	
        return diags([(b+c/s)/s*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])	
        

    return f,g,h

def f_logDate():
    def f(x,*args):
        return sum([log(abs(y))**2 for y in x[:-1]])

    def g(x,*args):
        return np.array([2*log(abs(z))/z for z in x[:-1]] + [0])

    def h(x,*args):
        #return np.diag([(2-2*log(abs(y)))/y**2 for y in x[:-1]]+[0])	
        return diags([(2-2*log(abs(y)))/y**2 for y in x[:-1]]+[0])	

    return f,g,h

def f_logDate_sqrt_b(pseudoCount=1e-8):
    def f(x,*args):
        return sum([sqrt(b)*log(abs(y))**2 if b>0 else pseudoCount*log(abs(y))**2 for (y,b) in zip(x[:-1],args[0])])

    def g(x,*args):
        return np.array([2*sqrt(b)*log(abs(z))/z if b>0 else 2*pseudoCount*log(abs(z))/z  for (z,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        #return np.diag([sqrt(b)*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])	
        return diags([sqrt(b)*(2-2*log(abs(y)))/y**2 if b>0 else pseudoCount*(2-2*log(abs(y)))/y**2  for (y,b) in zip(x[:-1],args[0])]+[0])	

    return f,g,h


def f_logDate_log_b():
    def f(x,*args):
        return sum([log(1+sqrt(b))*log(abs(y))**2 for (y,b) in zip(x[:-1],args[0])])

    def g(x,*args):
        return np.array([2*log(1+sqrt(b))*log(abs(z))/z for (z,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        #return np.diag([log(1+sqrt(b))*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])
        return diags([log(1+sqrt(b))*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])

    return f,g,h    
    
def logIt(tree,smpl_times,root_age=None,seqLen=1000,brScale=None,c=10,x0=None,f_obj=None,maxIter=MAX_ITER,min_b=0):
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2
    #N = len([node for node in tree.postorder_node_iter() if node is not tree.seed_node and node.edge_length > min_b])
    cons_mtrx = []
    cons_residue = []
    
    idx = 0
    b = [1.]*N

    for node in tree.postorder_node_iter():
        #if node is not tree.seed_node and node.edge_length > min_b:
        node.idx = idx
        idx += 1
        #else:
        #    node.idx = -1    
        if node.is_leaf():
            node.height = 0
            node.constraint = [0.0]*(N+1)
            node.residue = 0
            node.constraint[node.idx] = node.edge_length
            if node.edge_length > min_b:
                #node.constraint[node.idx] = node.edge_length
                b[node.idx] = node.edge_length
            else:
                #node.residue = node.edge_length
                b[node.idx] = 0
            node.constraint[N] = -smpl_times[node.taxon.label]
        else:
            children = list(node.child_node_iter())           
            node.height = min(children[0].height,children[1].height)
            a = [ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N+1) ]
            r = children[1].residue - children[0].residue
            cons_mtrx.append(a)
            cons_residue.append(r)

            if node is not tree.seed_node: 
                if children[0].height < children[1].height:
                    node.constraint = children[0].constraint
                    node.residue = children[0].residue
                else:
                    node.constraint = children[1].constraint
                    node.residue = children[1].residue
                node.constraint[node.idx] = node.edge_length
                if node.edge_length > min_b:
                    #node.constraint[node.idx] = node.edge_length
                    b[node.idx] = node.edge_length
                else:
                    #node.residue += node.edge_length
                    b[node.idx] = 0
            elif root_age is not None:
                a = children[0].constraint[:-1] + [children[0].constraint[-1]+root_age]
                r = -children[0].residue
                cons_mtrx.append(a)    
                cons_residue.append(r)

    x0 = ([1.]*N + [0.01]) if x0 is None else x0
    bounds = Bounds(np.array([MIN_NU]*N + [MIN_RATE]),np.array([np.inf]*(N+1)),keep_feasible=True)
    args = (b)
    linear_constraint = LinearConstraint(csr_matrix(cons_mtrx),[0]*len(cons_mtrx),[0]*len(cons_mtrx))
    #linear_constraint = LinearConstraint(csr_matrix(cons_mtrx),cons_residue,cons_residue)
    #linear_constraint = LinearConstraint(cons_mtrx,[0]*len(cons_mtrx),[0]*len(cons_mtrx))


    if f_obj is not None:
        f,g,h = f_obj(x,args)
    elif brScale == 'sqrt':
        f,g,h = f_logDate_sqrt_b(pseudoCount=min_b)
    elif brScale == 'log':
        f,g,h = f_logDate_log_b()
    elif brScale == 'lsd':
        f,g,h = f_logDate_lsd(c=c,s=seqLen)    
    else:
        f,g,h = f_logDate()    

    result = minimize(fun=f,method="trust-constr",x0=x0,bounds=bounds,args=args,constraints=[linear_constraint],options={'disp': True,'verbose':3,'maxiter':maxIter},jac=g,hess=h)
    x = result.x
    mu = x[N]
    fx = f(x,args)
        
    return mu,fx,x


def logDate_with_random_init(tree,sampling_time=None,root_age=None,leaf_age=None,brScale=False,nrep=1,min_nleaf=3,seqLen=1000,maxIter=MAX_ITER,seed=None,min_b="AUTO"):
    if min_b is None:
        min_b = 0
    elif min_b == "AUTO":
        min_b = minVar_bisect([node.edge_length for node in tree.postorder_node_iter() if node is not tree.seed_node])     
    
    print("Pseudo: " + str(min_b))

    smpl_times = {}
    
    if sampling_time is None:
        if leaf_age is None:
            leaf_age = 1
        for node in tree.leaf_node_iter():
            smpl_times[node.taxon.label] = leaf_age
        if root_age is None:
            root_age = 0    
    else:        
        with open(sampling_time,"r") as fin:
            fin.readline()
            for line in fin:
                name,time = line.split()
                smpl_times[name] = float(time)
    X,seed = random_date_init(tree,smpl_times,nrep,min_nleaf=min_nleaf,rootAge=root_age,seed=seed,min_b=0)
    print("Finished initialization with random seed " + str(seed))
    f_min = None
    x_best = None

    for i,x0 in enumerate(X):
        _,f,x = logIt(tree,smpl_times,root_age=root_age,brScale=brScale,x0=x0,seqLen=seqLen,maxIter=maxIter,min_b=min_b)
        if f_min is None or f < f_min:
            f_min = f
            x_best = x
            s_tree,t_tree = scale_tree(tree,x_best)
            print("Found a better log-scored configuration")
            print("New mutation rate: " + str(x_best[-1]))
            print("New log score: " + str(f_min))
            print("Scaled tree")
            print(s_tree.as_string("newick"))
            print("Time tree")
            print(t_tree.as_string("newick"))
             
    
    return x_best[-1],f_min,x_best,s_tree,t_tree 
    

def logDate_with_lsd(tree,sampling_time,root_age=None,brScale=False,lsdDir=None,seqLen=1000,maxIter=MAX_ITER,min_b="AUTO"):
    if min_b is None:
        min_b = 0
    elif min_b == "AUTO":
        min_b = minVar_bisect([node.edge_length for node in tree.postorder_node_iter() if node is not tree.seed_node])     
    
    wdir = run_lsd(tree,sampling_time,outputDir=lsdDir)
    
    x0 = read_lsd_results(wdir)
    x1 = [1.]*len(x0)
    #x1 = [1.]*(len([node for node in tree.postorder_node_iter() if node is not tree.seed_node and node.edge_length >= min_b])+1)
    x1[-1] = x0[-1] 
    #x1=x0
    smpl_times = {}

    with open(sampling_time,"r") as fin:
        fin.readline()
        for line in fin:
            name,time = line.split()
            smpl_times[name] = float(time)


    #mu,f,x = calibrate_log_opt(tree,smpl_times,root_age=root_age,brScale=brScale,x0=x1)
    mu,f,x = logIt(tree,smpl_times,root_age=root_age,brScale=brScale,x0=x1,seqLen=seqLen,maxIter=maxIter,min_b=min_b)
    
    if lsdDir is None:
        rmtree(wdir)

    s_tree,t_tree = scale_tree(tree,x) 

    return mu,f,x,s_tree,t_tree   

def run_lsd(tree,sampling_time,outputDir=None):
    wdir = outputDir if outputDir is not None else mkdtemp()
    treefile = normpath(join(wdir,"mytree.tre"))
    print(treefile)
    tree_swift = treeswift.read_tree_dendropy(tree)
    #tree.write_to_path(treefile,"newick")
    tree_swift.write_tree_newick(treefile)
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
    cons_mtrx = []
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
            cons_mtrx.append({'type':'eq','fun':g,'args':(a,)})

            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length
                b[node.idx] = node.edge_length
            elif root_age is not None:
                a = np.array(children[0].constraint[:-1] + [children[0].constraint[-1]-root_age])
                cons_mtrx.append({'type':'eq','fun':g,'args':(a,),'jac':g_gradient})    

    x0 = ([1.]*N + [0.01]) if x0 is None else x0
    bounds = [(0.00000001,999999)]*(N+1)
    args = (b)
    
    if brScale:
        result = minimize(fun=f2,x0=x0,args=args,bounds=bounds,constraints=cons_mtrx,method="SLSQP",options={'maxiter': 1500, 'ftol': 1e-10,'disp': True,'iprint':3},jac=f2_gradient)
    else:
        result = minimize(fun=f1,x0=x0,bounds=bounds,constraints=cons_mtrx,method="SLSQP",jac=f1_gradient ,options={'maxiter': 1500, 'ftol': 1e-10,'disp': True,'iprint':3})
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
    mu = x[-1]

    for node in tree.postorder_node_iter():
        if node is not tree.seed_node and node.idx >= 0:
            key = node.bipartition    
            mapping[key] = node.idx

    count = 0
    for node in s_tree.postorder_node_iter():
        if node is not s_tree.seed_node:
            if node.bipartition in mapping:
                idx = mapping[node.bipartition]
                node.edge_length *= x[idx]
                count += 1

    t_tree = Tree.get(data=s_tree.as_string("newick"),taxon_namespace=taxa,schema="newick",rooting="force-rooted")
    
    for node in t_tree.postorder_node_iter():
        if node is not t_tree.seed_node:
            node.edge_length /= mu

    return s_tree,t_tree        
