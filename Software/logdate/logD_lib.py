from dendropy import Tree,TaxonNamespace
import numpy as np
#from QP import quadprog_solve_qp, cvxopt_solve_qp
from math import exp,log, sqrt
from scipy.stats.stats import pearsonr
from scipy.optimize import minimize
from os.path import basename, dirname, splitext,realpath,join,normpath,isdir,isfile,exists

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
    x0 = [0]*N + [mu]
    
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
            if el > 0:
                x0[idx] = node.edge_length/float(el)

    return x0        
        

def calibrate_log_opt(tree,smpl_times,root_age=None,brScale=False,x0=None):
    def f0(x,*args):
        return sum([b*(w-1)*(w-1) for (w,b) in zip(x[:-1],args[0])])

    def f1(x):
        return sum([log(y)*log(y) for y in x[:-1]])
    
    def f2(x,*args):
        return sum([b*log(y)*log(y) for (y,b) in zip(x[:-1],args[0])])

    def g(x,a):    
        return sum([ x[i]*a[i] for i in range(len(x)) ])

    def h(x,p):
        return x[p]

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
                cons_eq.append({'type':'eq','fun':g,'args':(a,)})    

    x0 = ([1.]*N + [0.01]) if x0 is None else x0
    bounds = [(0.00000001,999999)]*(N+1)
    args = (b)
    #x1 = minimize(fun=f0,x0=x0,args=args,bounds=bounds,constraints=cons_eq+cons_ineq,method="SLSQP").x
    
    #x0 = [1.]*N + [0.02]
    #bounds = [(0.00000001,999999)]*(N+1)
    
    if brScale:
        result = minimize(fun=f2,x0=x0,args=args,bounds=bounds,constraints=cons_eq,method="SLSQP")
    else:
        result = minimize(fun=f1,x0=x0,bounds=bounds,constraints=cons_eq,method="SLSQP")
    x = result.x
    s = x[N]
    
    
    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
            node.edge_length *= x[node.idx]/s
    
    '''
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            #node.time = node.smplTime
            node.time = -node.constraint[N]
            #print(node.taxon.label + " " + str(node.time))
        else:
            children = list(node.child_node_iter())
            node.time = children[0].time - x[children[0].idx]*children[0].edge_length/s
            #print(node.label +  " " + str(node.time))
        #if node is not tree.seed_node:
        #    node.edge_length *= x[node.idx]/s
    '''

    return s

