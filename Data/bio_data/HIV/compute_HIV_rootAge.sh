#! /bin/bash

tree=$1
smplTime=$2

rootTime=`compute_nodeAge.py -i $tree -t $smplTime -m root_age`
date -d @`echo "3600*24*365*$rootTime" | bc -l` 
