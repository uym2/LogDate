#! /usr/bin/env python

# suppose the trees have the same topology but different branch lengths

from dendropy import Tree, TaxonNamespace

def main():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--input",required=True,help="Input trees")
    parser.add_argument("-o","--output",required=True,help="Output trees")
    args = vars(parser.parse_args())

    taxa = TaxonNamespace()
    filein = args["input"]
    fileout = args["output"]
    acc_brlen = {}
    nTrees = 0
    # Although using TreeList provided in Dendropy can be a more convenient solution,
    # I opted out for that because it requires storing a large number of trees in the memory at the same time
    # If the input trees are big then we will run out of memory 
    # Had problem with a set of 7k trees of 10k leaves which required >60G of memory just to store the trees
    # Here I read each tree and process it one-by-one. 
    #Just have to be thoughtful about making the taxon_namespace shared among all the trees
    print("Reading...")
    with open(filein,'r') as fin:
        strings = fin.readlines()       
        for s in strings:
            tree = Tree.get(data=s,schema="newick",taxon_namespace=taxa,rooting="force-rooted")
            tree.encode_bipartitions()
            for node in tree.preorder_node_iter():
                if node is not tree.seed_node:
                    key = node.bipartition
                    if not key in acc_brlen:
                        acc_brlen[key] = node.edge_length
                    else:
                        acc_brlen[key] += node.edge_length
            nTrees += 1
    
    print("Computing new tree ...")
    summarized_tree = tree
    for node in summarized_tree.preorder_node_iter():
        key = node.bipartition
        if key in acc_brlen:
            node.edge_length = acc_brlen[key]/float(nTrees)
    
    print("Writing...")
    with open(fileout,'w') as fout:
        fout.write(summarized_tree.as_string("newick"))    
                                   
if __name__=="__main__":
    main()                                
