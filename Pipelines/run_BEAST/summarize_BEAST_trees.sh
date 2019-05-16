BEAST_trees=$1
nTrees=$2
sample_trees=$3 # output
summarized_tree=$4 # output

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

tail -n $((nTrees+1)) $BEAST_trees | head -n $nTrees  | sed -e "s/tree.*\] //g" | sed  "s/\[\&rate=[0-9\.]*\]//g" > $sample_trees

temp=`mktemp`

python $DIR/summarize_trees.py -i $sample_trees -o $temp

# add trailing zero: just for this data set. The trees from LSD have trailing zeros for leaf names but BEAST removed them. Here we add the trailing zeros back
python $DIR/add_trailing_zeros.py $temp $summarized_tree 6

rm $temp
