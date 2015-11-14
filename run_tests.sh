#!/bin/bash
#BSUB -W 20:10                # wall-clock time (hrs:mins)
#BSUB -e cluster-errors.%J    # error file name in which %J is replaced by the job ID
#BSUB -o clutser-output.%J    # output file name in which %J is replaced by the job ID
#BSUB -n 1

# Run R non-interactively.
R CMD BATCH run_tests.R
