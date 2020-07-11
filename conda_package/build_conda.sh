#! /bin/bash

conda-build logdate/
conda convert --platform all //anaconda3/conda-bld/osx-64/logdate-1.5.0*.tar.bz2 -o //anaconda3/conda-bld/
anaconda upload --force //anaconda3/conda-bld/*/logdate-1.5.0*.tar.bz2
