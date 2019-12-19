#! /usr/bin/env python

import treeswift

tree_lsd = treeswift.read_tree_newick("phyML_noOG_rooted.tre.result.date.newick")
tree_logD = treeswift.read_tree_newick("wlogDu001q10r.tre")
tree_phyML_Bstr = treeswift.read_tree_newick("phyML_B_strict.tre")
tree_phyML_Bln = treeswift.read_tree_newick("phyML_B_lnorm.tre")
tree_Bstr = treeswift.read_tree_newick("B_strict.tre")
tree_Bln = treeswift.read_tree_newick("B_lnorm.tre")
tree_lf = treeswift.read_tree_newick("LF_1.tre")


ltt_lsd = tree_lsd.ltt(present_day=2011.44, show_plot=False)
ltt_logD = tree_logD.ltt(present_day=2011.44, show_plot=False)
ltt_phyML_Bstr = tree_phyML_Bstr.ltt(present_day=2011.44, show_plot=False)
ltt_phyML_Bln = tree_phyML_Bln.ltt(present_day=2011.44, show_plot=False)
ltt_Bstr = tree_Bstr.ltt(present_day=2011.44, show_plot=False)
ltt_Bln = tree_Bln.ltt(present_day=2011.44, show_plot=False)
ltt_lf = tree_lf.ltt(present_day=2011.44, show_plot=False)

data = [(ltt_lsd,"LSD"),(ltt_lf,"LF"),(ltt_logD,"wLogDate"),(ltt_phyML_Bstr,"phyML_B_strict"),(ltt_phyML_Bln,"phyML_B_lnorm"),(ltt_Bstr,"B_strict"),(ltt_Bln,"B_lnorm")]


for D,name in data:
    for t in D:
        print(name + " " + str(t) + " " + str(D[t]))        
