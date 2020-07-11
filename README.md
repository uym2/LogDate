# Installation
If you have Anaconda, you can install LogDate with conda install

``` bash
   conda install -c uym2 logdate 
```  

Alternatively, you can either clone the repository to your machine 
```bash
   git clone https://github.com/uym2/LogDate.git
```
or simply download the zip file. To install, go to the LogDate folder and type

``` bash
  python setup.py develop
```

# Usage
```bash
launch_LogDate.py [-h] -i INPUT [-t SAMPLINGTIME] -o OUTPUT [-v]
                         [-p REP] [-s RSEED] [-l SEQLEN] [-m MAXITER]
                         [-u ADDPSEUDO] [-z ZERO]
```
Arguments include

```bash
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input trees
  -t SAMPLINGTIME, --samplingTime SAMPLINGTIME
                        Sampling times / Calibration points. Default: root = 0
                        and all leaves = 1
  -o OUTPUT, --output OUTPUT
                        Output trees
  -v, --verbose         Show verbose message. Default: NO
  -p REP, --rep REP     The number of random replicates for initialization.
                        Default: use 1 initial point
  -s RSEED, --rseed RSEED
                        Random seed to generate starting tree initial points
  -l SEQLEN, --seqLen SEQLEN
                        The length of the sequences. Default: 1000
  -m MAXITER, --maxIter MAXITER
                        The maximum number of iterations for optimization.
                        Default: 50000
  -u ADDPSEUDO, --addpseudo ADDPSEUDO
                        Add pseudo count for per-branch weighting. Default:
                        0.01
  -z ZERO, --zero ZERO  Set zero-length branches (if any) to this number.
                        LogDate cannot process zero-length branches. Default:
                        1e-10
```

# Examples
```bash
python launch_LogDate.py -i tests/D750_11_10_1_phyML_lnorm.tre -o myTimeTree.tre -t tests/D750_11_10_1_smplTime_1.txt 
```
