#!/bin/bash
#$ -N full
#$ -q LT
#$ -pe smp 4
#$ -wd /Users/ssrivastva/lswasp/ml/code/
#$ -t 1-10

module load R/3.3.2_gcc-5.4.0

R CMD BATCH --no-save --no-restore "--args 0 $SGE_TASK_ID" submit.R samp/full_ml_$SGE_TASK_ID.rout
