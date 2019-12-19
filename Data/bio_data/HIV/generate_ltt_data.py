#! /usr/bin/env python

import treeswift

tree_lsd = treeswift.read_tree_newick("lsd_ingroups.tre")
tree_logD = treeswift.read_tree_newick("wlogDr10u001_ingroups.tre")
tree_lf = treeswift.read_tree_newick("LF_ingroups.tre")


ltt_lsd = tree_lsd.ltt(present_day=2018.08, show_plot=False)
ltt_logD = tree_logD.ltt(present_day=2018.08, show_plot=False)
ltt_lf = tree_lf.ltt(present_day=2018.08, show_plot=False)

data = [(ltt_lsd,"LSD"),(ltt_lf,"LF"),(ltt_logD,"wLogDate")]


for D,name in data:
    for t in D:
        print(name + " " + str(t) + " " + str(D[t]))        
