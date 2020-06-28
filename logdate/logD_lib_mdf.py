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
import treeswift
import cvxpy as cp

MAX_ITER = 50000
MIN_NU = 1e-12
MIN_MU = 1e-5


#lsd_file = "../lsd-0.2/bin/lsd.exe" if platform.system() == "Linux" else "../lsd-0.2/src/lsd"
#lsd_exec=normpath(join(dirname(realpath(__file__)),lsd_file))

def f_logDate_sqrt_scale(pseudo=0,seqLen=1000):
# square-root scaling strategy: let nu'_i = sqrt(b_i)*nu_i; solve the optimzation problem for nu' instead of nu
    def f(x,*args):
        return sum([log(abs(y/sqrt(b)))**2 for (y,b) in zip(x[:-1],args[0])])

    def g(x,*args):
        return np.array([2*(log(abs(y/sqrt(b))))/y for (y,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        return diags([2*(1-log(abs(y/sqrt(b))))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])	

    return f,g,h

def f_wlogDate_sqrt_scale(pseudo=0,seqLen=1000):
# square-root scaling strategy: let nu'_i = sqrt(b_i)*nu_i; solve the optimzation problem for nu' instead of nu
    def f(x,*args):
        w = [ sqrt(b+pseudo/seqLen) for b in args[0]]
        return sum([a*(log(abs(y/sqrt(b))))**2 for (y,a,b) in zip(x[:-1],w,args[0])])

    def g(x,*args):
        w = [ sqrt(b+pseudo/seqLen) for b in args[0]]
        return np.array([a*2*(log(abs(y/sqrt(b))))/y for (y,a,b) in zip(x[:-1],w,args[0])] + [0])

    def h(x,*args):
        w = [ sqrt(b+pseudo/seqLen) for b in args[0]]
        return diags([a*2*(1-log(abs(y/sqrt(b))))/y**2 for (y,a,b) in zip(x[:-1],w,args[0])]+[0])	

    return f,g,h

def f_logDate_linear_scale(pseudo=0,seqLen=1000):
# square-root scaling strategy: let nu'_i = sqrt(b_i)*nu_i; solve the optimzation problem for nu' instead of nu
    def f(x,*args):
        return sum([log(abs(y/b))**2 for (y,b) in zip(x[:-1],args[0])])

    def g(x,*args):
        return np.array([2*(log(abs(y/b)))/y for (y,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        return diags([2*(1-log(abs(y/b)))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])	

    return f,g,h

def f_wlogDate_linear_scale(pseudo=0,seqLen=1000):
# square-root scaling strategy: let nu'_i = sqrt(b_i)*nu_i; solve the optimzation problem for nu' instead of nu
    def f(x,*args):
        w = [ sqrt(b+pseudo/seqLen) for b in args[0]]
        return sum([a*(log(abs(y/b)))**2 for (y,a,b) in zip(x[:-1],w,args[0])])

    def g(x,*args):
        w = [ sqrt(b+pseudo/seqLen) for b in args[0]]
        return np.array([a*2*(log(abs(y/b)))/y for (y,a,b) in zip(x[:-1],w,args[0])] + [0])

    def h(x,*args):
        w = [ sqrt(b+pseudo/seqLen) for b in args[0]]
        return diags([a*2*(1-log(abs(y/b)))/y**2 for (y,a,b) in zip(x[:-1],w,args[0])]+[0])	

    return f,g,h

def f_lsd(pseudo=10,seqLen=1000):
    def f(x,*args):
        return sum([(b+pseudo/seqLen)*(y-1)*(y-1) for (y,b) in zip(x[:-1],args[0])])

    def g(x,*args):
        return np.array([2*(b+pseudo/seqLen)*(y-1) for (y,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        return diags([2*(b+pseudo/seqLen) for (y,b) in zip(x[:-1],args[0])]+[0])	

    return f,g,h

def f_LF(pseudo=0,seqLen=1000):
    def f(x,*args):
        return sum([b*(y-log((abs(y)))) for (y,b) in zip(x[:-1],args[0])])

    def g(x,*args):
        return np.array([b*(1-1/y) for (y,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        return diags([b/(y*y) for (y,b) in zip(x[:-1],args[0])]+[0])	

    return f,g,h

def f_PL(alpha_lf=1,alpha_lg=1):
# PL = "penalized likelihood"
    def f(x,*args):
        return sum([alpha_lf*b*(y-log(abs(y))) + alpha_lg*sqrt(b)*log(abs(y))**2 for (y,b) in zip(x[:-1],args[0])])
    
    def g(x,*args):
        return np.array([alpha_lf*b*(1-1/abs(y)) + alpha_lg*sqrt(b)*2*log(abs(y))/abs(y) for (y,b) in zip(x[:-1],args[0])]+[0])
        
    def h(x,*args):
        return diags([alpha_lf*b/(y*y) + alpha_lg*sqrt(b)*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-1],args[0])] + [0])

    return f,g,h

def f_logDate(pseudo=0,seqLen=1000):
    def f(x,*args):
        return sum([log(abs(y))**2 for y in x[:-1]])

    def g(x,*args):
        return np.array([2*log(abs(z))/z for z in x[:-1]] + [0])

    def h(x,*args):
        #return np.diag([(2-2*log(abs(y)))/y**2 for y in x[:-1]]+[0])	
        return diags([(2-2*log(abs(y)))/y**2 for y in x[:-1]]+[0])	

    return f,g,h

def f_logDate_sqrt_b(pseudo=0,seqLen=1000):
    def f(x,*args):
        return sum([sqrt(b+pseudo/seqLen)*log(abs(y))**2 for (y,b) in zip(x[:-1],args[0])])

    def g(x,*args):
        return np.array([2*sqrt(b+pseudo/seqLen)*log(abs(z))/z for (z,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        return diags([sqrt(b+pseudo/seqLen)*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])	

    return f,g,h


def f_logDate_log_b(pseudo=0,seqLen=1000):
    def f(x,*args):
        return sum([log(1+sqrt(b+pseudo/seqLen))*log(abs(y))**2 for (y,b) in zip(x[:-1],args[0])])

    def g(x,*args):
        return np.array([2*log(1+sqrt(b+pseudo/seqLen))*log(abs(z))/z for (z,b) in zip(x[:-1],args[0])] + [0])

    def h(x,*args):
        return diags([log(1+sqrt(b+pseudo/seqLen))*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[:-1],args[0])]+[0])

    return f,g,h    


def setup_constraint(tree,smpl_times,root_age=None,scale=None):
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2
    cons_eq = []
    
    idx = 0
    b = [1.]*N
   
    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1
        new_constraint = None        
        lb = node.taxon.label if node.is_leaf() else node.label
        if lb in smpl_times:
            node.height = 0
            new_constraint = [0.0]*(N+1)
            new_constraint[N] = -smpl_times[lb]            
        
        if not node.is_leaf():
            c1,c2 = list(node.child_node_iter()) # assuming each internal node has exactly two children
            f0 = lb in smpl_times
            f1 = c1.constraint is not None
            f2 = c2.constraint is not None

            if f1 and f2:
                if not f0:
                    node.height = min(c1.height,c2.height) + 1
                    new_constraint = c1.constraint if c1.height < c2.height else c2.constraint
                    a = [ (c1.constraint[i] - c2.constraint[i]) for i in range(N+1) ]
                    cons_eq.append(a)
                else:    
                    a1 = c1.constraint[:-1] + [-smpl_times[lb]-c1.constraint[-1]]
                    a2 = c2.constraint[:-1] + [-smpl_times[lb]-c2.constraint[-1]] 
                    cons_eq.append(a1)
                    cons_eq.append(a2)
            elif f1:
                if not f0:
                    node.height = c1.height
                    new_constraint = c1.constraint
                else:    
                    a1 = c1.constraint[:-1] + [-smpl_times[lb]-c1.constraint[-1]]
                    cons_eq.append(a1)
            elif f2:
                if not f0:
                    node.height = c2.height
                    new_constraint = c2.constraint        
                else:        
                    a2 = c2.constraint[:-1] + [-smpl_times[lb]-c2.constraint[-1]] 
                    cons_eq.append(a2)   
            
        if new_constraint is not None:    
            node.constraint = new_constraint
        if node is not tree.seed_node:    
            node.constraint[node.idx] = node.edge_length
            b[node.idx] = node.edge_length
    
    return cons_eq,b
                        

def setup_constraint_old(tree,smpl_times,root_age=None,scale=None):
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2
    cons_eq = []
    
    idx = 0
    b = [1.]*N

    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1
        if node.is_leaf():
            node.height = 0
            node.constraint = [0.0]*(N+1)
            if scale is None:
                node.constraint[node.idx] = node.edge_length
            elif scale == 'sqrt':    
                node.constraint[node.idx] = sqrt(node.edge_length) 
            elif scale == 'linear':
                node.constraint[node.idx] = 1 
            node.constraint[N] = -smpl_times[node.taxon.label]
            b[node.idx] = node.edge_length
        else:
            children = list(node.child_node_iter())           
            node.height = min(children[0].height,children[1].height)+1
            a = [ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N+1) ]
            cons_eq.append(a)

            if node is not tree.seed_node: 
                node.constraint = children[0].constraint if children[0].height < children[1].height else children[1].constraint
                if scale is None:
                    node.constraint[node.idx] = node.edge_length
                elif scale == 'sqrt':    
                    node.constraint[node.idx] = sqrt(node.edge_length) 
                elif scale == 'linear':
                    node.constraint[node.idx] = 1 
                b[node.idx] = node.edge_length
            elif root_age is not None:
                a = children[0].constraint[:-1] + [children[0].constraint[-1]-root_age]
                cons_eq.append(a)    

    return cons_eq,b

    
def logIt(tree,smpl_times,f_obj,scale=None,root_age=None,x0=None,maxIter=MAX_ITER,pseudo=0,seqLen=1000):
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2

    cons_eq,b = setup_constraint(tree,smpl_times,root_age=root_age,scale=scale)

    if scale is None:
        bounds = Bounds(np.array([MIN_NU]*N+[MIN_MU]),np.array([np.inf]*(N+1)),keep_feasible=False)
        x_init = [1.]*N + [0.01] if x0 is None else x0
    elif scale == 'sqrt':
        bounds = Bounds(np.array([MIN_NU*sqrt(a) for a in b]+[MIN_MU]),np.array([np.inf]*(N+1)),keep_feasible=False)
        x_init = [sqrt(a) for a in b] + [0.01] if x0 is None  else [sqrt(a)*x for (x,a) in zip (x0[:-1],b)]+[x0[-1]]
    elif scale == 'linear':    
        bounds = Bounds(np.array([MIN_NU*a for a in b]+[MIN_MU]),np.array([np.inf]*(N+1)),keep_feasible=False)
        x_init = [a]*N + [0.01] if x0 is None else [a*x for (x,a) in zip (x0[:-1],b)]+[x0[-1]]
    
    args = (b)
    linear_constraint = LinearConstraint(csr_matrix(cons_eq),[0]*len(cons_eq),[0]*len(cons_eq),keep_feasible=False)

    f,g,h = f_obj(pseudo=pseudo,seqLen=seqLen)
    
    print("Initial state:" )
    print("mu = " + str(x_init[-1]))
    print("fx = " + str(f(x_init,args)))
    print("Maximum constraint violation: " + str(np.max(csr_matrix(cons_eq).dot(x_init))))
    
    result = minimize(fun=f,method="trust-constr",x0=x_init,bounds=bounds,args=args,constraints=[linear_constraint],options={'disp': True,'verbose':3,'maxiter':maxIter},jac=g,hess=h)
   
    if scale is None:
        x_opt = result.x
    elif scale == 'sqrt':
         x_opt = [ x/sqrt(a) for (x,a) in zip(result.x[:-1],b) ] + [result.x[-1]]
    elif scale == 'linear':     
         x_opt = [ x/a for (x,a) in zip(result.x[:-1],b) ] + [result.x[-1]]

    mu = x_opt[N]
    fx = result.fun  
    
    return mu,fx,x_opt  

def setup_smpl_time(tree,sampling_time=None,root_age=None,leaf_age=None):
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
    
    return smpl_times   
    
def find_LF_opt(cons_eq,b):
    C = [ a[:-1] for a in cons_eq]
    #C = [ [ (-1)**(b<0)*sqrt(abs(b)) for b in a[:-1] ] for a in cons_eq]
    dt = [a[-1] for a in cons_eq]
    N = len(b)

    x = cp.Variable(N) # nu
    y = cp.Variable()  # mu

    param = cp.Parameter((1,len(b)),nonneg=True,value=np.reshape(np.array(b),(1,len(b))))
    #objective = cp.Minimize(param*(x-cp.log(x)))
    objective = cp.Minimize(cp.sqrt(param)*x-param*cp.log(x))
    constraints = [np.array([MIN_NU]*N) <= x, MIN_MU <= y, csr_matrix(C)*x+y*np.array(dt)==0]
    prob = cp.Problem(objective, constraints)
    fx = prob.solve(verbose=True) #,mosek_params={'MSK_DPAR_BASIS_REL_TOL_S':1e-20,'MSK_DPAR_BASIS_TOL_S':1e-9,'MSK_DPAR_BASIS_TOL_X':1e-9})
    x_opt = np.append([ x/sqrt(a) for (x,a) in zip (x.value,b) ],y.value)
    
    return x_opt,fx

def run_LF_cvxpy(tree,sampling_time=None,root_age=None,leaf_age=None):
    smpl_times = setup_smpl_time(tree,sampling_time=sampling_time,root_age=root_age,leaf_age=leaf_age)
    
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2

    cons_eq,b = setup_constraint(tree,smpl_times,root_age=root_age,sqrt_scale=True)    
    x_opt,fx = find_LF_opt(cons_eq,b)    
    s_tree,t_tree = scale_tree(tree,x_opt)
             
    return x_opt[-1],fx,x_opt,s_tree,t_tree 

def logDate_with_penalize_llh(tree,sampling_time=None,root_age=None,leaf_age=None,maxIter=MAX_ITER):
    smpl_times = setup_smpl_time(tree,sampling_time=sampling_time,root_age=root_age,leaf_age=leaf_age)
    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2

    cons_eq,b = setup_constraint(tree,smpl_times,root_age=root_age,sqrt_scale=True)    
    x_lf,_ = find_LF_opt(cons_eq,b)   
   
    '''
    max_nu = 10 #max(x_lf[:-1])
    print("max_nu = " + str(max_nu))
    
    alpha1 = log(max_nu)*log(max_nu)
    alpha2 = (max_nu - log(max_nu)) - 1
    print("alpha_lf = " + str(alpha1) + ", alpha_lg = " + str(alpha2))
    '''
    alpha1 = alpha2 = 1

    _,f_opt,x_opt = logIt(tree,smpl_times,f_PL,sqrt_scale=False,root_age=root_age,x0=x_lf,maxIter=maxIter)
    s_tree,t_tree = scale_tree(tree,x_opt)
   
    return x_opt[-1],f_opt,x_opt,s_tree,t_tree
    
    

def logDate_with_random_init(tree,f_obj,sampling_time=None,root_age=None,leaf_age=None,nrep=1,min_nleaf=3,maxIter=MAX_ITER,seed=None,scale=None,pseudo=0,seqLen=1000):
    smpl_times = setup_smpl_time(tree,sampling_time=sampling_time,root_age=root_age,leaf_age=leaf_age)
    
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            node.fixed_age = smpl_times[node.taxon.label]
        else:    
            node.fixed_age = None
    
    X,seed,_ = random_date_init(tree,smpl_times,nrep,min_nleaf=min_nleaf,rootAge=root_age,seed=seed)
    print("Finished initialization with random seed " + str(seed))
    f_min = None
    x_best = None

    i = 0
    n_succeed = 0

    for i,x0 in enumerate(X):
        x0 = X[i]
        _,f,x = logIt(tree,smpl_times,f_obj,root_age=root_age,x0=x0,maxIter=maxIter,scale=scale,pseudo=pseudo,seqLen=seqLen)
        print("Found local optimal for Initial point " + str(i+1))
        n_succeed += 1                
        
        if f_min is None or f < f_min:
            f_min = f
            x_best = x
            s_tree,t_tree = scale_tree(tree,x_best)
            print("Found a better log-scored configuration")
            print("New mutation rate: " + str(x_best[-1]))
            print("New log score: " + str(f_min))
            #print("Scaled tree")
            #print(s_tree.as_string("newick"))
            print("Time tree")
            print(t_tree.as_string("newick"))
             
    return x_best[-1],f_min,x_best,s_tree,t_tree 
    

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
    mu = x[-1]

    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
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
