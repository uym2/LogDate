from bitsets import bitset
from dendropy import Tree,Node
from sys import argv
from math import log, sqrt, exp
import random

EPSILON_nu = 1e-5
EPSILON_t = 1e-3

def random_date_init(tree, sampling_time, rep, rootAge=None, min_nleaf=3, seed=None):
    if seed is None:
        seed = random.randint(0,1024)

    random.seed(a=seed)    
    initial_fix_age(tree)
    
    history = []
    X = []
    T0 = []
    
    node_list = get_node_list(tree,min_nleaf=min_nleaf)
    k = len(node_list)  
    node_bitset = bitset('node_bitset',node_list)

    for i in range(rep):
        x = int(random.getrandbits(k))

        while x <= 0 or x in history:
            x = int(random.getrandbits(k))

        history.append(x)
        selected = list(node_bitset.fromint(x))
        x0,t0 = date_as_selected(tree,sampling_time,selected,rootAge=rootAge)
        X.append(x0)
        T0.append(t0)
    
    return X,seed,T0   
         
def print_date(tree,fout=None):
    for node in tree.postorder_node_iter():
       if node is not tree.seed_node:
           label = node.taxon.label if node.is_leaf() else node.label
           if fout is not None:
               fout.write(label + " " + str(node.time - node.parent_node.time) + "\n")
           else:
               print(label + " " + str(node.time - node.parent_node.time))
                     

def get_node_list(tree, min_nleaf=3):
    node_list = []
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.nleaf = 1
        else:
            node.nleaf = sum([x.nleaf for x in node.child_node_iter()])
            if node is not tree.seed_node and node.nleaf >= min_nleaf and node.fixed_age is None:
                node_list.append(node)
    return tuple(node_list)         

def date_by_RTT(tree,sampling_time,rootAge=None,epsilon_t=EPSILON_t):
    # preprocessing
    for node in tree.postorder_node_iter():
        node.as_leaf = False
        if node.is_leaf():
            node.time = sampling_time[node.taxon.label]
        else:
            node.time = None
    t0 = rootAge if rootAge is not None else compute_date_as_root(tree.seed_node)          
    
    if t0 is None:
        t_min = min(node.time for node in tree.preorder_node_iter() if node.time is not None)
        t0 = t_min - epsilon_t
    
    preprocess_node(tree.seed_node)

    tree.seed_node.time = t0
    date_from_root_and_leaves(tree.seed_node)
        
    x = compute_mu_from_dated_tree(tree)     
    return x,t0


def date_as_selected(tree,sampling_time,selected,rootAge=None,epsilon_t=EPSILON_t):
# selected is a list of nodes; they MUST be sorted in postorder
# each node in tree must have fixed_age attribute and at least two of them MUST NOT be None
    # refresh initial values
    preprocess_tree(tree,sampling_time)

    t0 = rootAge if rootAge is not None else compute_date_as_root(tree.seed_node)

    # dating
    for node in selected:
        # date the clade that is rooted at node
        preprocess_node(node)
        t = compute_date_as_root(node,t0 = t0)
        if t is not None:
            node.time = t
            date_from_root_and_leaves(node)
            node.as_leaf = True
    
    t_min = min(node.time for node in tree.preorder_node_iter() if node.time is not None)
    t0 = t_min - epsilon_t
    
    preprocess_node(tree.seed_node)

    tree.seed_node.time = t0
    date_from_root_and_leaves(tree.seed_node)
        
    x = compute_mu_from_dated_tree(tree)     
    return x,t0

def preprocess_tree(tree,sampling_time):
# all nodes in tree must have 'fixed_age' attribute and at least two of them MUST NOT be None
# the sampling_time is not needed in this updated version, but is kept to be consistent with 
# the (obsolete) downstream codes that use this function   
   for node in tree.postorder_node_iter():
       node.time = node.fixed_age
       node.as_leaf = (node is not tree.seed_node and node.fixed_age is not None and node.parent_node.fixed_age is None)
       
def initial_fix_age(tree):
# heuristically compute the time of the nodes below fixed-age nodes    
# and fix their ages to this heuristic (not a great idea...)    
   for node in tree.postorder_node_iter():
       node.time = node.fixed_age
       if node.time is not None:
           preprocess_node(node)
           date_from_root_and_leaves(node)     
           node.as_leaf = True 
    
   for node in tree.postorder_node_iter():
       node.fixed_age = node.time

def preprocess_node(a_node):
    #print("Current node: " + a_node.label)
    stack = [(a_node,False)]

    while stack:
        node,passed = stack.pop()
        if node.is_leaf() or node.as_leaf:
           node.nearest_leaf = node
           node.nearest_t = node.time
           node.delta_b = node.edge_length
           lb = node.taxon.label if node.is_leaf() else node.label
           #print(lb)
        elif passed:
           min_t = None
           min_child = None
           for c in node.child_node_iter():
                if min_t is None or c.nearest_t < min_t or (c.nearest_t == min_t and c.delta_b < min_child.delta_b): 
                    min_t = c.nearest_t
                    min_child = c
           node.nearest_leaf = min_child.nearest_leaf
           node.nearest_t = min_child.nearest_t
           node.delta_b = min_child.delta_b + ( node.edge_length if node.edge_length else 0)
        else:
           stack.append((node,True))
           stack += [(x,False) for x in node.child_node_iter()]


def compute_date_as_root(a_node,t0=None):
    SD = 0
    ST = 0
    SDT = 0
    STS = 0
    n = 0
    t_min = None

    a_node.d2root = 0
    stack = a_node.child_nodes()
    while stack:
        node = stack.pop()
        node.d2root = node.parent_node.d2root + node.edge_length
        if node.is_leaf() or node.as_leaf:
            SDT += node.d2root*node.time
            SD  += node.d2root
            STS += node.time**2
            ST  += node.time
            n   += 1
            t_min = node.time if t_min is None else min(t_min,node.time)
        else:
            stack += node.child_nodes()

    t = None
    if SD*ST != n*SDT:
        t_star = (STS*SD - SDT*ST) / (SD*ST - n*SDT)
        if t_star < t_min and (t0 is None or t_star > t0):
            t = t_star
    return t
                
def date_from_root_and_leaves(root_node):
    stack = root_node.child_nodes()
    while stack:
        node = stack.pop()
        if not node.is_leaf() and not node.as_leaf:
            stack = stack + node.child_nodes()
        if node.time is not None:
            assert node.time > root_node.time
            continue
        u = node.nearest_leaf
        delta_t = node.nearest_t - node.parent_node.time
        mu = node.delta_b / delta_t

        while u != node:
            v = u.parent_node
            v.time = u.time - u.edge_length/mu
            assert v.time < u.time
            assert v.time > root_node.time
            u = v  
        assert node.time > root_node.time    

            
            
def compute_mu_from_dated_tree(tree):
    a = 0
    b = 0
    c = 0
    x0 = []
    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
            if (node.time == node.parent_node.time):
                continue
            alpha = node.edge_length / (node.time - node.parent_node.time)
            w = log(1 + sqrt(node.edge_length))
            a += w
            b -= 2*w*log(alpha)
            c += w*(log(alpha)**2)

    mu = exp(-b/2/a)
  
    for node in tree.postorder_node_iter():
        if node is not tree.seed_node:
            nu = mu*(node.time - node.parent_node.time)/node.edge_length if (node.time != node.parent_node.time) else EPSILON_nu
            x0.append(nu)

    x0.append(mu)        
            
    return x0      

def main():
   tree_file = argv[1]
   sampling_time_file = argv[2]

   tree = Tree.get_from_path(tree_file,"newick")
   sampling_time = {}

   with open(sampling_time_file,'r') as fin:
       fin.readline()
       for line in fin:
           taxon, time = line.split()
           sampling_time[taxon] = float(time)
   
   random_date_init(tree, sampling_time, 10, min_nleaf=8)

   '''        
   for node in tree.postorder_node_iter():
       node.as_leaf = False
       if node.is_leaf():
           node.time = sampling_time[node.taxon.label]
           node.nearest_leaf = node
           node.nearest_t = node.time
           node.delta_b = node.edge_length
       else:
           node.time = None
           min_t = None
           min_child = None
           for c in node.child_node_iter():
                if min_t is None or c.nearest_t < min_t or (c.nearest_t == min_t and c.delta_b < min_child.delta_b): 
                    min_t = c.nearest_t
                    min_child = c
           node.nearest_leaf = min_child.nearest_leaf
           node.nearest_t = min_child.nearest_t
           node.delta_b = min_child.delta_b + node.edge_length
            


   compute_date_as_root(tree.seed_node)
   date_from_root_and_leaves(tree.seed_node)
   
   for node in tree.postorder_node_iter():
       if node is not tree.seed_node:
           label = node.taxon.label if node.is_leaf() else node.label
           print(label + " " + str(node.time - node.parent_node.time))
   '''
   #print("mu = " + str(compute_mu_from_dated_tree(tree))) 
if __name__ =="__main__":
    main()                     
