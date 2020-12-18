#!/bin/bash
#$ -N real_samp500_100
#$ -q LT
#$ -pe smp 8
#$ -l h_rt=520:00:00
#$ -l s_rt=520:00:00
#$ -wd /Users/ssrivastva/real_dls/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-1000
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/
#$ -j y

module load R/3.3.2_gcc-5.4.0

R CMD BATCH --no-save --no-restore "--args 11 $SGE_TASK_ID" submit.R samp/real_samp500_k100$SGE_TASK_ID.rout

