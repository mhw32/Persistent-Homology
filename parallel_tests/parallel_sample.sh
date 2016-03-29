#!/bin/bash
#BSUB -W 01:00                # wall-clock time (hrs:mins)
#BSUB -e error.err   # error file name in which %J is replaced by the job ID
#BSUB -o output.out    # output file name in which %J is replaced by the job ID
#BSUB -q shared
#BSUB -J parallel_sample
#BSUB -u mike.wu@yale.edu

module load Apps/R
module load Rpkgs/RMPI
module load Rpkgs/DOPARALLEL

cd /home/fas/cisewski/mhw32/scratch/homology/parallel_tests
mpirun -n 1 R --slave -f parallel_sample.R
