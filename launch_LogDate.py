import logdate
from logdate.logD_lib import logDate_with_random_init,f_wLogDate
from dendropy import Tree
from os import remove,path
from logdate.tree_lib import tree_as_newick
import argparse
from sys import argv

print("Launching " + logdate.PROGRAM_NAME + " version " + logdate.PROGRAM_VERSION)
print(logdate.PROGRAM_NAME + " was called as follow")
print(" ".join(argv))

parser = argparse.ArgumentParser()

parser.add_argument("-i","--input",required=True,help="Input trees")
parser.add_argument("-t","--samplingTime",required=False,help="Sampling times / Calibration points. Default: root = 0 and all leaves = 1")
parser.add_argument("-b","--backward",action='store_true',help="Use backward time. Default: NO")
parser.add_argument("-o","--output",required=True,help="Output trees")
parser.add_argument("-v","--verbose",action='store_true',help="Show verbose message. Default: NO")
parser.add_argument("-k","--keeplabel",action='store_true',help="Suppress auto label assignment to internal nodes. WARNING: LogDate uses the labels to identify calibration points. Use this option only if your tree has UNIQUE LABELING FOR ALL nodes. Default: NO")
parser.add_argument("-p","--rep",required=False,help="The number of random replicates for initialization. Default: use 1 initial point")
parser.add_argument("-s","--rseed",required=False,help="Random seed to generate starting tree initial points")
parser.add_argument("-l","--seqLen",required=False,help="The length of the sequences. Default: 1000")
parser.add_argument("-m","--maxIter",required=False,help="The maximum number of iterations for optimization. Default: 50000")
parser.add_argument("-u","--addpseudo",required=False,help="Add pseudo count for per-branch weighting. Default: 0.01")
parser.add_argument("-z","--zero",required=False,help="Set zero-length branches (if any) to this number. LogDate cannot process zero-length branches. Default: 1e-10")

args = vars(parser.parse_args())

sampling_time = args["samplingTime"]
nrep = int(args["rep"]) if args["rep"] else 1
seqLen = int(args["seqLen"]) if args["seqLen"] else 1000
pseudo = 0.01 if args["addpseudo"] is None else float(args["addpseudo"])
maxIter = int(args["maxIter"]) if args["maxIter"] else 50000
randseed = int(args["rseed"]) if args["rseed"] else None
zero_len = float(args["zero"]) if args["zero"] else 1e-10

f_obj = f_wLogDate
verbose = args["verbose"]
do_label = not args["keeplabel"]
bw_time = args["backward"]

with open(args["input"],'r') as fin:
    tree_strings = fin.readlines()

if path.exists(args["output"]):
    remove(args["output"])

for treestr in tree_strings:  
    tree = Tree.get(data=treestr,schema='newick',preserve_underscores=True)
    # labeling
    if do_label:
        nodeIdx = 0
        for node in tree.preorder_node_iter():
            if not node.is_leaf():
                node.label = "I" + str(nodeIdx)
                nodeIdx += 1
    # handle zero-length branches
    for node in tree.preorder_node_iter():            
        if node is not tree.seed_node and node.edge_length == 0:
            node.edge_length = zero_len
    # dating        
    mu,f,x,s_tree,t_tree = logDate_with_random_init(tree,f_obj,sampling_time,bw_time=bw_time,nrep=nrep,min_nleaf=10,maxIter=maxIter,seed=randseed,pseudo=pseudo,seqLen=seqLen,verbose=verbose)
    tree_as_newick(t_tree,outfile=args["output"],append=True)
    print("Clock rate: " + str(mu))
    print("Log score: " + str(f))
