#!/bin/bash
#$ -N ml
#$ -q INFORMATICS
#$ -pe smp 4
#$ -wd /Users/ssrivastva/lswasp/ml/code/
#$ -t 1-1000

module load R/3.3.2_gcc-5.4.0

R CMD BATCH --no-save --no-restore "--args 1 $SGE_TASK_ID" submit.R samp/ml_$SGE_TASK_ID.rout
