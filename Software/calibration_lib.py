from dendropy import Tree
import numpy as np
#from QP import quadprog_solve_qp, cvxopt_solve_qp
from math import exp,log, sqrt
from scipy.stats.stats import pearsonr
from scipy.optimize import minimize

def lnorm_deviation_from_clock(tree,sd,rate):
# Note: we force the distribution to have mean 1, so 
# there is only 1 parameter to control the lognormal distribution
# sd here is the standard deviation of the lognormal distribution, 
# NOT its underlying normal distribution

    tree1 = Tree(tree)
    mu = -0.5*log(sd*sd+1)
    sigma = sqrt(log(sd*sd+1))

    for node in tree1.postorder_node_iter():
        if node is tree1.seed_node:
            continue
        f = np.random.lognormal(mean=mu,sigma=sigma)
        node.edge_length = node.edge_length*f*rate
    return tree1    

def exp_deviation_from_clock(tree,rate):
# Note: we force the distribution to have mean 1, so 
#there is no free parameter for an exponential distribution

    tree1 = Tree(tree)

    for node in tree1.postorder_node_iter():
        if node is tree1.seed_node:
            continue
        f = np.random.exponential()
        node.edge_length = node.edge_length*f*rate
    return tree1    

def gamma_deviation_from_clock(tree,shape,rate):
    tree1 = Tree(tree)

    for node in tree1.postorder_node_iter():
        if node is tree1.seed_node:
            continue
        f = np.random.gamma(shape,scale=1/shape)
        node.edge_length = node.edge_length*f*rate
    return tree1    


def deviation_from_clock(tree):
    factor = []
    for node in tree.postorder_node_iter():
        if node is tree.seed_node:
            continue
        f = np.random.gamma(exp(1.5),scale=1/exp(1.5))
        #f = np.random.lognormal(1,0.4)
        factor.append(f)
        node.edge_length = node.edge_length*f

    return factor    

def calibrate_composite_opt(tree,smpl_times,root_age=None,c=10.0,L=1219.0):
    def f0(x):
        return sum([log(y)*log(y) for y in x[:-1]])

    def f1(x,*args):
        v = args[1]
        return sum([log(w)*log(w)/v + L*b*(w-1)*(w-1)/w for (w,b) in zip(x[:-1],args[0])])
        #return sum([L*b*(w-1)*(w-1)/w for (w,b) in zip(x[:-1],args[0])])
        #return sum([b*b*(w-1)*(w-1)/(b+c/L) for (w,b) in zip(x[:-1],args[0])])
        #return sum([1.0/y*1.0/y-2*1.0/y for y in x[:-1]])

    def g(x,a):    
        return sum([ x[i]*a[i] for i in range(len(x)) ])


    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2
    cons = []
    
    idx = 0
    b = [1.]*N

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
            cons.append({'type':'eq','fun':g,'args':(a,)})

            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length
                b[node.idx] = node.edge_length
            elif root_age is not None:
                a = np.array(children[0].constraint[:-1] + [children[0].constraint[-1]-root_age])
                cons.append({'type':'eq','fun':g,'args':(a,)})    

    x0 = [1.]*N + [0.01]
    bounds = [(0.00000001,999999)]*(N+1)
    x1 = minimize(fun=f0,x0=x0,bounds=bounds,constraints=cons,method="SLSQP").x

    l_sum = 0
    l_ssq = 0

    for y in x1[:-1]:
        l_sum += log(y)
        l_ssq += log(y)*log(y)

    v = l_ssq/N - (l_sum/N)*(l_sum/N)    
    print(v)
    
    b1 = [w_i*b_i for (w_i,b_i) in zip(x1[:-1],b)]

    args = (b1,v)
    bounds = [(0.00000001,999999)]*(N+1)
    
    x2 = minimize(fun=f1,x0=x1,args=args,bounds=bounds,constraints=cons,method="SLSQP").x
    
    s = x2[N]
    
    
    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
            node.edge_length *= x2[node.idx]/s
    
    return s

def calibrate_log_opt(tree,smpl_times,root_age=None,brScale=False):
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

    x0 = [1.]*N + [0.01]
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

def calibrate_with_sampling_time(tree,smpl_times,verbose=False):
    if verbose:
        print("Start calibrating tree. Setting up the matrices and constraints...")

    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2


    h = np.array([0.0]*(N+1)).reshape((N+1,))
    b = np.array([0.0]*(n-1)).reshape((n-1,)) 
    #b = np.array([0.0]*(n-1)).reshape((n-1,)) 
    #b = []

    G = np.negative(np.identity(N+1))
    #P = N*np.identity(N) - np.ones((N,N))
    #q = np.array([0.0]*N).reshape((N,))

    P = np.identity(N+1)
    P[N][N] = 0
    q = np.array([-2.0]*N+[0]).reshape((N+1,)) 
    
    #A = [1.0]*N
    A = np.zeros((n-1,N+1))
    #A[0] = [1.0]*N + [0]
    #A = None

    idx = 0
    r = 0
    c = 0
    L = 1000.0

    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1

         
        if node is not tree.seed_node: 
            P[node.idx][node.idx] = node.edge_length*node.edge_length/(node.edge_length+c/L)
            q[node.idx] = -2.0*node.edge_length*node.edge_length/(node.edge_length+c/L)
        

        if node.is_leaf():
            node.constraint = [0.0]*(N+1)
            node.constraint[node.idx] = node.edge_length
            node.constraint[N] = -smpl_times[node.taxon.label]
            #node.smplTime = smpl_times[node.taxon.label]
        else:
            children = list(node.child_node_iter())
                       
            a = [ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N+1) ]

            #A = np.vstack([A,a])
            #A = np.vstack([A,a]) if A is not None else a
            #b.append(children[0].smplTime - children[1].smplTime)
            A[r] = a
            r += 1

            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length
                #node.smplTime = children[0].smplTime

    if verbose:
        print("Done with setting up. Started quadratic optimization ...")

    f = cvxopt_solve_qp(P,q,G,h,A,b)
    #f = quadprog_solve_qp(P,q,G,h,A,b)
    w = f[N]
    print("Clock rate: " + str(w))
    if verbose:
        print("Optimal solution found. Calibrating tree ...")
    
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            #node.time = node.smplTime
            node.time = -node.constraint[N]
            #print(node.taxon.label + " " + str(node.time))
        else:
            children = list(node.child_node_iter())
            node.time = children[0].time - f[children[0].idx]*children[0].edge_length/w
            #print(node.label +  " " + str(node.time))
        #if node is not tree.seed_node:
        #    node.edge_length *= f[node.idx]/s
    
    return f


def calibrate_tree(tree,verbose=False):
    if verbose:
        print("Start calibrating tree. Setting up the matrices and constraints...")

    n = len(list(tree.leaf_node_iter()))
    N = 2*n-2

    h = np.array([0.0]*N).reshape((N,))
    #b = np.array([N]+[0.0]*(n-1)).reshape((n,)) 
    b = np.array([0.0]*(n-1)).reshape((n-1,)) 

    G = np.negative(np.identity(N))
    P = np.identity(N)
    #q = np.array([0.0]*N).reshape((N,))
    q = np.array([-2.0]*N).reshape((N,)) 
    
    #A = [1.0]*N
    #A = np.zeros((n,N))
    #A[0] = [1.0]*N
    A = None

    idx = 0
    r = 1
    for node in tree.postorder_node_iter():
        node.idx = idx
        idx += 1
        if node.is_leaf():
            node.constraint = [0.0]*N
            node.constraint[node.idx] = node.edge_length
        else:
            children = list(node.child_node_iter())
                       
            a = [ (children[0].constraint[i] - children[1].constraint[i]) for i in range(N) ]
            #A = np.vstack([A,a])
            A = np.vstack([A,a]) if A is not None else a
            #A[r] = a
            r += 1
            if node is not tree.seed_node: 
                node.constraint = children[0].constraint
                node.constraint[node.idx] = node.edge_length

    if verbose:
        print("Done with setting up. Started quadratic optimization ...")

    #f = cvxopt_solve_qp(P,q,G,h,A,b)
    f = quadprog_solve_qp(P,q,G,h,A,b)
    if verbose:
        print("Optimal solution found. Calibrating tree ...")

    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
            node.edge_length *= f[node.idx]

    return f

def main():
    from sys import argv

    tree = Tree.get_from_path(argv[1],"newick")
    
    '''
    smpl_times = {}
    with open(argv[2],"r") as fin:
        fin.readline()
        for line in fin:
            name,time = line.split()
            smpl_times[name] = float(time)
    '''
    f = deviation_from_clock(tree)

    tree.write_to_path(argv[2],"newick")
   
    '''
    m=sum([1/x for x in f])/len(f)
    with open('f.txt','w') as fout:
        for x in f:
            fout.write(str(1/x/m) + "\n")
    '''    

    #f = calibrate_with_sampling_time(tree,smpl_times)
    #f = calibrate_tree(tree)
    #print(f)
    '''
    m = sum(f1)/len(f1)
    with open('f1.txt','w') as fout:
        for x in f1:
            fout.write(str(x/m) + "\n")
    '''

    #tree.write_to_path(argv[2],"newick")
    
    #print(np.corrcoef([1/x for x in f],f1))

if __name__=="__main__":
    main()     
