from dendropy import Tree
import numpy as np
from scipy.sparse import csr_matrix,diags
from scipy.optimize import LinearConstraint, minimize,Bounds
from logdate.util_lib import minVar_bisect
from logdate.init_lib import random_date_init, preprocess_tree, date_from_root_and_leaves, preprocess_node, date_by_RTT
from math import sqrt,log
from sys import stdout

MAX_ITER = 50000
MIN_NU = 1e-12
MIN_MU = 1e-5
EPSILON_t = 1e-5

def write_time_tree(tree,outfile=None,append=False):
    if outfile:
        outstream = open(outfile,'a' if append else 'w')
    else:
        outstream = stdout

    __write__(tree.seed_node, outstream)
    outstream.write(";\n")
    if outfile:
        outstream.close()

def __write__(node, outstream):
    if node.is_leaf():
        try:
            outstream.write(node.taxon.label)
        except:
            outstream.write(node.label)
    else:
        outstream.write('(')
        is_first_child = True
        for child in node.child_node_iter():
            if is_first_child:
                is_first_child = False
            else:
                outstream.write(',')
            __write__(child,outstream)
        outstream.write(')')
        if node.label is not None:
            outstream.write(str(node.label))
    
    if not node.parent_node is None:
        outstream.write(":" + str(node.time-node.parent_node.time))

def f_wlogDate():
    def f(x,*args):
        return sum([sqrt(b)*log(abs(y))**2 for (y,b) in zip(x[2:],args[0])])

    def g(x,*args):
        return np.array([0,0]+[2*sqrt(b)*log(abs(z))/z for (z,b) in zip(x[2:],args[0])])

    def h(x,*args):
        return diags([0,0]+[sqrt(b)*(2-2*log(abs(y)))/y**2 for (y,b) in zip(x[2:],args[0])])	

    return f,g,h

def logging(f_obj,constraints,weights,x0,maxIter=MAX_ITER):
    f,g,h = f_obj()
    N = len(x0)
    bounds = Bounds(np.array([-np.inf,MIN_MU]+[MIN_NU]*(N-2)),np.array([np.inf]*N))
    result = minimize(fun=f,method="trust-constr",x0=x0,bounds=bounds,args=(weights,),constraints=[constraints,],options={'disp': True,'verbose':3,'maxiter':maxIter},jac=g,hess=h)
    return result

def log_from_random_init(tree,sampling_time,root_age=None,leaf_age=None,brScale=False,nrep=1,min_nleaf=3,maxIter=MAX_ITER,seed=None,max_cutoff=1e-3,seqLen=1000):
    # mark short terminal branches
    x,_ = date_by_RTT(tree,sampling_time,rootAge=root_age)
    mu_RTT = x[-1]  
    tauMin = EPSILON_t
    tauMax = max(tauMin,0.1/seqLen/mu_RTT)
    print("Short-branch range estimation: ")     
    print("mu_RTT = " + str(mu_RTT))
    print("tauMin =  " + str(tauMin) + " tauMax = " + str(tauMax))
    L = [node.edge_length for node in tree.postorder_node_iter() if node is not tree.seed_node]
    cutoff = min(max_cutoff,minVar_bisect(L))
    calibs,count_short = calibs_from_leaf_times(tree,sampling_time,short_terms_thres=cutoff,tauMin=tauMin,tauMax=tauMax)
    constrs,weights = setup_constraints(tree,calibs)
    print("Finished setting up constraints according to sampling time. Cutoff threshold for short branches set to " + str(cutoff) + ". Eliminated " + str(count_short) + " short terminal branches and relaxed their constraints.")

    # multiple random initial points
    X,seed,T0 = random_date_init(tree,sampling_time,nrep,min_nleaf=min_nleaf,rootAge=root_age,seed=seed)
    print("Finished initialization of " + str(nrep) + " replicates with random seed " + str(seed))
    
    f_min = None
    x_best = None

    # hacking (dirty) code to match the 'init_idx' and 'constr_idx'
    idx = 0
    for node in tree.postorder_node_iter():
        node.init_idx = idx
        idx += 1 
    
    for i,x1 in enumerate(X):
        print("Start searching from initial point " + str(i+1)) 
        
        # matching constr_idx and init_idx
        x0 = [0]*(len(weights)+2)
        for node in tree.postorder_node_iter():
            if node.constr_idx is not None:
                x0[node.constr_idx] = x1[node.init_idx]
        x0[0] = T0[i]*x1[-1]       # mu*t0
        mu = x0[1] = x1[-1]        # mu
        
        print("Initial state:")
        print("mu = " + str(x1[-1]))
        print("t0 = " + str(T0[i]))
                
        result = logging(f_wlogDate,constrs,weights,x0,maxIter=maxIter)
        x = result.x
        f = result.fun

        if f_min is None or f < f_min:
            f_min = f
            x_best = x
            compute_time_tree(tree,x_best,sampling_time)
            print("Found a better log-scored solution")
            print("New mutation rate: " + str(x_best[1]))
            print("New root age: " + str(x_best[0]/x_best[1]))
            print("New log score: " + str(f_min))
            print("New time tree")
            write_time_tree(tree)
    return x_best,f_min

def compute_time_tree(tree,x_best,sampling_time):
    preprocess_tree(tree,sampling_time)
    preprocess_node(tree.seed_node)

    mu = x_best[1]
    t0 = x_best[0]/mu
    tree.seed_node.time = t0

    unprocess_clades = []

    for node in tree.preorder_node_iter():
        if node is not tree.seed_node and node.constr_idx is not None:
            nu = x_best[node.constr_idx]
            node.time = node.parent_node.time + node.edge_length*nu/mu
            if node.has_short_child and not node.is_short:
                unprocess_clades.append(node)
    
    for node in unprocess_clades:            
        preprocess_node(node)
        date_from_root_and_leaves(node)

    #for node in tree.preorder_node_iter():
    #    if node is not tree.seed_node:
    #        node.edge_length = node.time - node.parent_node.time
                

def mark_short_terminals(tree,threshold):
    count_short = 0
    for node in tree.postorder_node_iter():
        if node is tree.seed_node:
            continue
        if node.is_leaf():
            node.is_short = node.edge_length <= threshold                
            node.has_short_child = False
        else:
            node.is_short = node.edge_length <= threshold and ( sum(not c.is_short for c in node.child_node_iter()) == 0 )   
            node.has_short_child = (sum(c.is_short for c in node.child_node_iter()) > 0)
        count_short += int(node.is_short)
    return count_short        
        
def calibs_from_leaf_times(tree,sampling_time,short_terms_thres=0,tauMin=EPSILON_t,tauMax=100*EPSILON_t):
    count_short = mark_short_terminals(tree,short_terms_thres)

    calibs = []

    for node in tree.postorder_node_iter():
        node.fixed_age = None
        if node is tree.seed_node:
            continue
        if node.is_leaf():
            node.fixed_age = node.tmin = node.tmax = sampling_time[node.taxon.label]
            if not node.is_short:
               calibs.append((node,node.tmin,node.tmax)) 
            node.h =  0   
        else:
            node.tmax = min(c.tmax for c in node.child_node_iter())
            node.tmin = None
            node.h = max(c.h for c in node.child_node_iter() if c.is_short) + 1 if node.has_short_child else 0
            if node.has_short_child and not node.is_short:
                print("Found one short clade with h = " + str(node.h))
                node.tmax -= node.h*tauMin
                node.tmin = node.tmax - node.h*tauMax
                calibs.append((node,node.tmin,node.tmax))      
                node.fixed_age = (node.tmin + node.tmax)/2                  
    return calibs,count_short    

def setup_constraints(tree,calibs):
    idx = 2
    row = 0
    
    for node in tree.postorder_node_iter():
        node.constr_idx = None

    rows = []
    cols = []
    data = []
    uppers = []
    lowers = []
    weights = []

    for cnode,tmin,tmax in calibs:
        if (tmin is None and tmax is None) or (tmin is not None and tmax is not None and tmin > tmax):
            print("Warning: invalid calibration found for tmin = " + str(tmin) + " and tmax = " + str(tmax))
            continue

        node = cnode
        double_constr_flag = (tmin is not None) and (tmax is not None) and (tmin < tmax)
        upper = 0 if tmax is not None else np.inf
        lower = 0 if tmin is not None else -np.inf

        # mu*t0
        rows.append(row)
        cols.append(0)
        data.append(1)

        # mu
        rows.append(row)
        cols.append(1)
        data.append(-tmin if tmin is not None else -tmax)
        
        if tmin is None:
            lowers.append(-np.inf)
            uppers.append(0)
        else:
            lowers.append(lower)
            uppers.append(upper if not double_constr_flag else np.inf)    
        
        if double_constr_flag:
            # mu*t0
            rows.append(row+1)
            cols.append(0)
            data.append(1)

            # mu
            rows.append(row+1)
            cols.append(1)
            data.append(-tmax)

            lowers.append(-np.inf)
            uppers.append(upper)
        
        while node is not tree.seed_node:
            if node.constr_idx is None:
                node.constr_idx = idx
                weights.append(node.edge_length)
                idx += 1        
            rows.append(row)
            cols.append(node.constr_idx)
            data.append(node.edge_length)
            if double_constr_flag:
                rows.append(row+1)
                cols.append(node.constr_idx)
                data.append(node.edge_length)
            
            node = node.parent_node    
        
        row += (1 + int(double_constr_flag))

    A = csr_matrix((data,(rows,cols)),shape=(max(rows)+1,max(cols)+1))
    return LinearConstraint(A,lowers,uppers),weights

def main():
    from sys import argv

    tree = Tree.get_from_path(argv[1],'newick')
    sampling_time = {}

    with open(argv[2],'r') as fin:
        fin.readline()
        for line in fin:
            taxon,time = line.split()
            sampling_time[taxon] = float(time)

    x_best = log_from_random_init(tree,sampling_time)


if __name__ == "__main__":
    main()    
