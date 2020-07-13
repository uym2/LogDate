from dendropy import Tree

def find_LCAs(myTree,myQueries):
# find LCA for each of the list of nodes in myQueries
# return the label of the LCA node of each query
# based on the algorithm described here https://cp-algorithms.com/graph/lca.html
    def euler_tour():
    # perform the euler_tour (DFS) on the tree
    # output: 
    #   + E: the euler tour
    #   + F: the index in E where each node first occurs
    #   + H: the height (branch distance to root) of each node in E
        E = []
        H = {}
        F = {}
        def __traverse__(node,idx,h):
            lb = node.taxon.label if node.is_leaf() else node.label
            E.append(node)
            H[node] = h
            F[lb] = idx
            next_idx = idx+1
            for c in node.child_node_iter():
                next_idx = __traverse__(c,next_idx,h+1)
                E.append(node)
                next_idx += 1
            return next_idx    
        __traverse__(myTree.seed_node,0,1)
        return E,F,H
    
    def min_segment_tree(E,H):
    # build a min segment-tree
    # to query the minimum of any range of H
    # in O(logn)
        n = len(E)
        t = [None]*(4*n)
        def __build__(node,b,e):
            if b==e:
                t[node] = E[b]
            else:
                mid = (b+e)//2
                __build__(2*node,b,mid) 
                __build__(2*node+1,mid+1,e)
                l = t[2*node]
                r = t[2*node+1]
                t[node] = l if H[l] < H[r] else r
        __build__(1,0,n-1)     
        return t                          

    def query_segment_tree(t,q,E,F,H):
        L = min(F[a] for a in q)
        R = max(F[a] for a in q)
        def __query__(node,b,e,L,R):
            if (R < b or L > e):
                return None
            if (b >= L and e <= R):
                return t[node]    
            mid = (b+e)//2
            left = __query__(node*2,b,mid,L,R)
            right = __query__(node*2+1,mid+1,e,L,R)
            if left is None:
                return right
            if right is None:
                return left
            return left if H[left] < H[right] else right        

        return __query__(1,0,len(E)-1,L,R)

    E,F,H = euler_tour()
    t = min_segment_tree(E,H)
    myLCAs = []
    for q in myQueries:
        lca = query_segment_tree(t,q,E,F,H)
        myLCAs.append(lca)
    return myLCAs    

'''
from sys import argv
myTree = Tree.get_from_path(argv[1],"newick")
queries = []
with open(argv[2],'r') as fin:
    for line in fin:
        q,t = line.strip().split()
        queries.append(q.strip().split('=')[-1].split('+'))
LCAs = find_LCAs(myTree,queries)
nodeIdx = 0
for q,lca in zip(queries,LCAs):
    lb = lca.label if not lca.is_leaf() else lca.taxon.label
    if not lb:
        lca.label = "I" + str(nodeIdx)
        nodeIdx += 1
    if lca.is_leaf():
        print(q,lca.taxon.label)
    else:
        print(q,lca.label)  
'''        
