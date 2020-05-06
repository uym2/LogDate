from bitsets import bitset
from dendropy import Tree,Node
from sys import argv
from math import log, sqrt, exp
import random

EPSILON_nu = 1e-5
EPSILON_t = 1e-3

def init_calibrate(tree,sampling_time):
    for node in tree.postorder_node_iter():
        node.time = None
        lb = node.taxon.label if node.is_leaf() else node.label
        if lb in sampling_time:
            node.time = sampling_time[lb]


def mark_active(tree,sampling_time):
# a node is active if it has a sampling time or one of its descendants is active
    for node in tree.postorder_node_iter():
        node.is_active = False
        lb = node.taxon.label if node.is_leaf() else node.label
        if lb in sampling_time:
            node.is_active = True
        elif not node.is_leaf():
            for c in node.child_node_iter():
                if c.is_active:
                    node.is_active = True
                    break

def mark_as_leaf(tree):
# a node is marked as_leaf if it is active but none of its children is active
# must call mark_active before this function
    for node in tree.postorder_node_iter():
        if node.is_active:
            node.as_leaf = True
            for c in node.child_node_iter():
                if c.is_active:
                    node.as_leaf = False
                    break
        else:
            node.as_leaf = False  

def compute_tmin(tree):
# compute tmin for each node according to the sampling_time
# this is to set the constraint on the randomly selected nodes to be dated using RTT
# if RTT result t does not agree with the constraint tmin < t then 
# we will postpone the RTT date of that node
# if a node is calibrated (i.e. fixed) at time t, then we set tmin t for that node
# MUST be called AFTER reset and mark_active
    for node in tree.preorder_node_iter():
        if node.time is not None: # this node is calibrated
            node.tmin = node.time
        else:
            node.tmin = None if node is tree.seed_node else node.parent_node.tmin

def reset(tree,sampling_time):
    init_calibrate(tree,sampling_time)
    mark_as_leaf(tree)

def preprocess(tree,sampling_time):
    mark_active(tree,sampling_time)
    init_calibrate(tree,sampling_time)
    mark_as_leaf(tree)
    compute_tmin(tree)

def get_calibrated_nodes(tree,sampling_time):
# return the list of the calibrated internal nodes
    node_list = []
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            lb = node.label if not node.is_leaf() else node.taxon.label
            if lb in sampling_time:
                node_list.append(node)
    return node_list            
        
def get_uncalibrated_nodes(tree,sampling_time,min_nleaf=3):
# return the list of all uncalibrated nodes with at least min_nleaf active descended leaves
# NOTE: MUST call mark_active and mark_as_leaf before this function
    node_list = []
    for node in tree.postorder_node_iter():
        if node.as_leaf:
            # a leaf can never be in this list because it must either be calibrated or inactive
            node.n_active_leaf = 1 if node.is_active else 0
        else:
            node.n_active_leaf = sum([x.n_active_leaf for x in node.child_node_iter()])
            if node is not tree.seed_node and node.n_active_leaf >= min_nleaf and node.label not in sampling_time:
                node_list.append(node)
    return tuple(node_list)         

def calibrate_set(tree,node_list):
    for node in node_list:
        node.time = "awaited"
    for node in tree.postorder_node_iter():    
        if node.time is not None:
            preprocess_node(node)
        if node.time == "awaited":             
            t = compute_date_as_root(node)
            if t is not None:
                node.time = t # successfully calibrated, otherwise, just skip  
            else:
                node.time = None
        if node.time is not None and not node.as_leaf:
            date_from_root_and_leaves(node)
            node.as_leaf = True

    # if the root has not been calibrated, now we have to calibrate it
    if tree.seed_node.time is None:
        preprocess_node(tree.seed_node)
        t = compute_date_as_root(tree.seed_node)
        tree.seed_node.time = t if t is not None else tree.seed_node.tmax - EPSILON_t
        date_from_root_and_leaves(tree.seed_node)
    return tree.seed_node.time    

def preprocess_node(a_node):
    stack = [(a_node,False)]

    while stack:
        node,passed = stack.pop()
        if node.as_leaf:
           node.nearest_calibrated = node
           node.nearest_t = node.time
           node.tmax = node.time
           node.delta_b = node.edge_length
        elif passed:
           min_t = None
           min_child = None
           for c in node.child_node_iter():
               if not c.is_active:
                   continue
               if min_t is None or c.nearest_t < min_t or (c.nearest_t == min_t and c.delta_b < min_child.delta_b): 
                    min_t = c.nearest_t
                    min_child = c
           node.nearest_calibrated = min_child.nearest_calibrated
           node.nearest_t = min_t
           node.tmax = min_t
           node.delta_b = min_child.delta_b + ( node.edge_length if node.edge_length else 0)
        else:
           stack.append((node,True))
           stack += [(x,False) for x in node.child_node_iter() if x.is_active]
                
def compute_date_as_root(a_node):
# NOTE: a_node MUST have tmin and tmax attributes
# the return value t must satisfy tmin < t < tmax, otherwise, return None
    SD = 0
    ST = 0
    SDT = 0
    STS = 0
    n = 0

    a_node.d2root = 0
    stack = [ c for c in a_node.child_node_iter() if c.is_active ]
    while stack:
        node = stack.pop()
        node.d2root = node.parent_node.d2root + node.edge_length
        if node.as_leaf:
            SDT += node.d2root*node.time
            SD  += node.d2root
            STS += node.time**2
            ST  += node.time
            n   += 1
        else:
            stack += [ c for c in node.child_node_iter() if c.is_active ] 

    t = None
    if SD*ST != n*SDT:
        t_star = (STS*SD - SDT*ST) / (SD*ST - n*SDT)
        if (a_node.tmax is None or t_star < a_node.tmax) and (a_node.tmin is None or t_star > a_node.tmin):
            t = t_star
    return t

def date_from_root_and_leaves(root_node):
# compute time for all active nodes below the specified root_node, heuristically
# assumming the specified root_node has been calibrated (i.e. root_node.time is not None) and
# the root_node has some active nodes below it
# also assumming that the preprocess_node function has been applied to root_node, so that
# all active nodes below it already have the needed attributes (e.g. nearest_t, delta_b, delta_t, etc.) 
    stack = [c for c in root_node.child_node_iter() if c.is_active]
    while stack:
        node = stack.pop()
        if not node.as_leaf:
            stack += [c for c in node.child_node_iter() if c.is_active]
        if node.time is not None:
            assert node.time > root_node.time
            continue
        u = node.nearest_calibrated
        delta_t = node.nearest_t - node.parent_node.time
        mu = node.delta_b / delta_t

        while u != node:
            v = u.parent_node
            v.time = u.time - u.edge_length/mu
            #print("Calibrated node " + v.label + " using " + root_node.label)
            assert v.time < u.time
            assert v.time > root_node.time
            u = v  
        assert node.time > root_node.time    

def get_init_from_dated_tree(tree):
    a = 0
    b = 0
    c = 0
    x0 = []
    for node in tree.postorder_node_iter():
        if node is not tree.seed_node and node.is_active: # the root does not have edge_length; inactive nodes do not have time attributes
            alpha = node.edge_length / (node.time - node.parent_node.time)
            w = log(1 + sqrt(node.edge_length))
            a += w
            b -= 2*w*log(alpha)
            c += w*(log(alpha)**2)

    mu = exp(-b/2/a)
  
    for node in tree.postorder_node_iter():
        if node is not tree.seed_node and node.is_active:
            #nu = mu*(node.time - node.parent_node.time)/node.edge_length if node.is_active else 1
            nu = mu*(node.time - node.parent_node.time)/node.edge_length 
            x0.append(nu)
    x0.append(mu)        
    return x0      

def get_random_initial(tree,random_subset,sampling_time):
    #reset(tree,sampling_time)
    t0 = calibrate_set(tree,random_subset)
    #t0 = scale_from_calibrated(tree)
    #print_date(tree)
    x = get_init_from_dated_tree(tree)
    return x,t0

def print_date(tree):
    for node in tree.postorder_node_iter():
       if node.is_active:
           label = node.taxon.label if node.is_leaf() else node.label
           #print(label + " " + str(node.time))

def random_date_init(tree, sampling_time, rep, min_nleaf=3, seed=None):
    if seed is None:
        seed = random.randint(0,1024)

    random.seed(a=seed)    
    
    history = []
    X = []
    T0 = []
    
    preprocess(tree,sampling_time)
    node_list = get_uncalibrated_nodes(tree,sampling_time,min_nleaf=min_nleaf)
    if not node_list:
        print("Warning: few calibration points were given. Set min_nleaf to 1 to maximize the number of possible initials")
        node_list = get_uncalibrated_nodes(tree,sampling_time,min_nleaf=1)
    k = len(node_list) # k should always be larger than 0, by construction

    nrep = rep
    is_lacking = False
    if k < rep and 2**k <= rep: # include k < rep to avoid computing 2**k if k is obviously large
        nrep = 2**k
        is_lacking = True
        print("Warning: do not have enough calibration points/sampling time to construct " + str(rep) + " distinct initials. Reduced to " + str(nrep))
    node_bitset = bitset('node_bitset',node_list)

    for i in range(nrep):
        if is_lacking:
            x = i
        else:    
            x = int(random.getrandbits(k))
            while x <= 0 or x in history:
                x = int(random.getrandbits(k))

        history.append(x)
        random_subset = list(node_bitset.fromint(x))
        x0,t0 = get_random_initial(tree,random_subset,sampling_time)
        X.append(x0)
        T0.append(t0)
        reset(tree,sampling_time)
    
    return X,seed,T0   
