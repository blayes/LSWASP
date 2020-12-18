#!/bin/bash
#$ -N comp_comb_samp250
#$ -q LT
#$ -pe smp 4
#$ -l h_rt=520:00:00
#$ -l s_rt=520:00:00
#$ -wd /Users/ssrivastva/dls/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/
#$ -j y

module load R

R CMD BATCH --no-save --no-restore "--args 6 $SGE_TASK_ID" submit_samp250.R samp/comp_comb_samp250$SGE_TASK_ID.rout
