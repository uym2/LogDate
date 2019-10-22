#! /bin/bash

# a wrapper to run r8s's implementation of LF
# only accept 1 input tree; require a file that specifies that sampling time (at all leaves)
# following the structure of what lsd requires, i.e.
# ntaxa
# taxon_1 time_1
# taxon_2 time_2
# ...
# taxon_n time_n

in_tree=$1
sampling_time=$2
nsites=$3
in_r8s=$4 # Despite its name, this is an output of this script: the Nexus file that will be fed to r8s
out_r8s=$5 # The output of r8s
out_tree=$6 # The output (i.e. calibrated) tree, which should be the last line of $out_r8s if r8s finished successfully

echo "Preparing input for r8s ..."

echo "     Copying the input tree ..."
echo "#nexus" > $in_r8s
echo "begin trees;" >> $in_r8s
echo -n "tree myTree = " >> $in_r8s
cat $in_tree >> $in_r8s
echo "end;" >> $in_r8s

echo "     Generating r8s commands ..."
echo "begin r8s;" >> $in_r8s
echo "blformat lengths=persite nsites=$nsites ultrametric=no;" >> $in_r8s

echo "          Converting sampling time of each leaf to the corresponding 'Fixed Age' in r8s..."
maxTime=`tail -n+2 $sampling_time | awk '{print $2}' | numlist -max`

is_first_line=true

while read line; do
    if [ "$is_first_line" = true ]; then
        is_first_line=false
    else    
        echo -n "fixage taxon=" >> $in_r8s
        echo -n $line | awk '{printf $1;}' >> $in_r8s
        echo -n " age=" >> $in_r8s
        echo -n $line | awk -v m=$maxTime '{printf m-$2;}' >> $in_r8s
        echo ";" >> $in_r8s
    fi    
done < $sampling_time 

echo "collapse;" >> $in_r8s
echo "divtime method=LF;" >> $in_r8s
echo "describe plot=chrono_description;" >> $in_r8s;
echo "end;" >> $in_r8s

echo "Running r8s"
r8s -b -f $in_r8s > $out_r8s
tail -n1 $out_r8s | sed "s/tree myTree = //g" > $out_tree 
