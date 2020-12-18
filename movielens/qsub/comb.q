#!/bin/bash
#$ -N comb
#$ -q LT
#$ -pe smp 4
#$ -l h_rt=520:00:00
#$ -l s_rt=520:00:00
#$ -wd /Users/ssrivastva/lswasp/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/
#$ -j y

module load R/3.3.2_gcc-5.4.0

R CMD BATCH --no-save --no-restore "--args 2 $SGE_TASK_ID" submit.R samp/comb_$SGE_TASK_ID.rout


