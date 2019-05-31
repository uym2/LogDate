#! /usr/bin/env python

# Suppose the tree's branch lengths was estimated by a ML tool such as RAxML or PhyML with true topology
# Since RAxML (and perhaps phyML as well) removes the root even if the tree's topology was given
# We want to root the estimated tree based on the true rooted tree
# We suppose the true rooted tree (i.e. the reference tree) is bifurcating

from dendropy import Tree, TaxonNamespace
from dendropy.datamodel.treemodel import Bipartition
from sys import argv

def reroot_at_edge(tree, edge, length1, length2):
# the method provided by dendropy DOESN'T seem to work ...
    if not edge:
	return
    head = edge.head_node
    tail = edge.tail_node
    if not tail:
	return

    if (length2 == 0) and head.is_leaf():
	return

    new_root = tree.node_factory()
    
    tail.remove_child(head)

    new_root.add_child(head)
    head.parent_node = new_root
    head.edge_length=length2

    p = tail.parent_node
    l = tail.edge_length

    new_root.add_child(tail)
    tail.parent_node = new_root
    tail.edge_length = length1

     	

    if (tail is tree.seed_node):
	head = new_root


    #while tail is not tree.seed_node:
    while p is not None:
	head = tail
	tail = p
        p = tail.parent_node
	
        l1 = tail.edge_length
	tail.remove_child(head)

	head.add_child(tail)
        tail.parent_node = head
	tail.edge_length=l
	l = l1
	
        if tail is tree.seed_node:
            break
	
    # out of while loop: tail IS now tree.seed_node
    if tail.num_child_nodes() == 1:
	# merge the 2 branches of the old root and adjust the branch length
	#sis = [child for child in tail.child_node_iter()][0]
	sis = tail.child_nodes()[0]
	l = sis.edge_length
	tail.remove_child(sis)    
	head.add_child(sis)
	sis.edge_length = l + tail.edge_length
	head.remove_child(tail)
	#tail.remove_child(head)



inTree=argv[1]
outTree=argv[2]
refTree=argv[3]


taxa = TaxonNamespace()
t = Tree.get_from_path(inTree,'newick',taxon_namespace=taxa)
rt = Tree.get_from_path(refTree,'newick',taxon_namespace=taxa)

c1,c2 = rt.seed_node.child_nodes()
e1 = c1.edge_length
e2 = c2.edge_length

a_leaf = list(rt.postorder_node_iter())[0]
print(a_leaf.taxon.label)
reroot_at_edge(t,a_leaf.edge,a_leaf.edge_length/2,a_leaf.edge_length/2)

t.is_rooted = True
rt.is_rooted = True

t.encode_bipartitions()
rt.encode_bipartitions()

for node in t.postorder_node_iter():
    if node.bipartition == c1.bipartition:
        break        

reroot_at_edge(t,node.edge,e1*node.edge_length/(e1+e2),e2*node.edge_length/(e1+e2))

t.write_to_path(outTree,'newick')
