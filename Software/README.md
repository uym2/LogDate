# Installation
## Download
If you have git, you can simply clone the LogDate repository to your machine `git clone https://github.com/uym2/LogDate.git`. Otherwise, you can download the zip file to your machine.
## Prerequisites
- Anaconda with Python 3.7: highly recommended
- Scipy 1.3.1: usually included with Anaconda
- bitsets: https://pypi.org/project/bitsets/

# Usage
```bash
launch_LogDate.py [-h] -i INPUT [-j OBJ] [-t SAMPLINGTIME] [-r ROOTAGE]
                         [-f LEAFAGE] -o OUTPUT [-v] [-p REP] [-s RSEED]
                         [-l SEQLEN] [-m MAXITER] [-u ADDPSEUDO] [-z ZERO]
```
Arguments include

```bash
 -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input tree
  -j OBJ, --obj OBJ     Objective function. Either LogDate or wLogDate.
                        Default: wLogDate
  -t SAMPLINGTIME, --samplingTime SAMPLINGTIME
                        Sampling time at leaf nodes. Default: None
  -r ROOTAGE, --rootAge ROOTAGE
                        Root age. Can be used with either -f or -t, but not
                        both. Default: None if -t is specified else 0
  -f LEAFAGE, --leafAge LEAFAGE
                        Leaf age. To be used with root age to infer relative
                        time. Will be overried by -t if -t is specified.
                        Default: None if -t is specified else 1.
  -o OUTPUT, --output OUTPUT
                        The output trees with branch lengths in time unit.
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
                        Add pseudo counting for per-branch weighting.Default:
                        0.01
  -z ZERO, --zero ZERO  Set zero-length branches (if any) to this number.
                        LogDate cannot process zero-length branches. Default:
                        1e-10
```
