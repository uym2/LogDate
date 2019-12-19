#! /usr/bin/env python

import treeswift

tree_lsd = treeswift.read_tree_newick("lsd.tre")
tree_logD = treeswift.read_tree_newick("wlogDu001q10r.tre")
tree_beast = treeswift.read_tree_newick("BEAST.tre")
tree_lf = treeswift.read_tree_newick("LF.tre")


ltt_lsd = tree_lsd.ltt(present_day=2015.81, show_plot=False)
ltt_logD = tree_logD.ltt(present_day=2015.81, show_plot=False)
ltt_beast = tree_beast.ltt(present_day=2015.81, show_plot=False)
ltt_lf = tree_lf.ltt(present_day=2015.81, show_plot=False)

data = [(ltt_lsd,"LSD"),(ltt_beast,"BEAST"),(ltt_lf,"LF"),(ltt_logD,"wLogDate")]


for D,name in data:
    for t in D:
        print(name + " " + str(t) + " " + str(D[t]))        
