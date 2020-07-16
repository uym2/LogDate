# Installation
If you have Anaconda, you can install LogDate with conda install

``` bash
   conda install -c uym2 logdate 
```  

Alternatively, you can either clone the repository to your machine 
```bash
   git clone https://github.com/uym2/LogDate.git
```
or simply download the zip file to your machine. To install, go to the LogDate folder and type

``` bash
  python setup.py develop
```

# Usage
LogDate accepts both calibration points (hard constraint on divergence time) of internal nodes and sampling times at leaf nodes, and a mixture of them. Here we give examples for the three most common use-cases:

## Use case 1: infer the unit ultrametric tree
If there is no calibrations given, LogDate will infer the unit (depth 1) ultrametric tree from an input tree.

``` bash
   python launch_LogDate.py -i <INPUT_TREE> -o <OUTPUT_TREE>
```

We give an example in folder `examples/unit_time_tree`, inside which you can find the sampled input tree `phylogram.nwk` and the sampled output tree `chronogram_unit.nwk`.

```bash
   cd examples/unit_time_tree
   python launch_LogDate -i phylogram.nwk -o wLogDate.nwk -s 1105
```

Here we use `-s` to specify a random seed. With the correct random seed given, the output file `wLogDate.nwk` should be an exact replicate of `chronogram_unit.nwk`.

## Use case 2: infer the time tree from phylodynamics data
A typical use-case in virus phylogeny is to infer the time tree from a phylogram inferred from sequences and their sampling times (i.e. calibration points given at leaf nodes). LogDate reads the calibration points / sampling times from an input file via the `-t` option.

```bash
   python launch_LogDate -i <INPUT_TREE> -o <OUTPUT_TREE> -t <SAMPLING_TIMES>
```

### 2.1. Complete sampling times
An example is given in the folder `examples/virus_all_samplingTime`, inside which you will find the sampled input tree (`phylogram.nwk`), sampled file to specify the sampling times (`sampling_time.txt`), and the sampled output tree of LogDate (`chronogram.nwk`). In this example, we give LogDate all the sampling times for all leaves (i.e. complete sampling times). The file `sampling_time.txt` has two columns, which are simply the species names and their corresponding sampling times.

```bash
   cd examples/virus_all_samplingTime
   python launch_LogDate -i phylogram.nwk -o wLogDate.nwk -t sampling.time.txt -s 1105
```

Again, `-s` is used to specify the random seed. Your output (`wLogDate.nwk`) should be an exact replicate of `chronogram.nwk`.

### 2.2. Partial (missing) sampling times
LogDate allows missing sampling times for a subset of the leaves, as long as there exists at least one pair of leaves with different sampling times. The usage of LogDate for is the same as in the case of complete sampling times. An example is given in the folder `example/viruse_some_samplingTime`. Here we specify the sampling times for only 52 species out of 110 total.

```bash
   cd examples/virus_some_samplingTime
   python launch_LogDate -i phylogram.nwk -o wLogDate.nwk -t sampling.time.txt -s 1105
```

### 2.3. Sampling times at internal nodes
LogDate allows the sampling times to be given both internal nodes and at leaves. An example is given in the folder `example/virus_internal_smplTime`. Here we use the node label to identify the internal nodes (see `example/virus_internal_smplTime/sampling_time.txt`). 
