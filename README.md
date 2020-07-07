# Installation
## Download
If you have git, you can simply clone the LogDate repository to your machine `git clone https://github.com/uym2/LogDate.git`. Otherwise, you can download the zip file to your machine.
## Prerequisites
- Anaconda with Python 3.7: highly recommended
- Scipy 1.3.1: usually included with Anaconda
- bitsets: https://pypi.org/project/bitsets/

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
