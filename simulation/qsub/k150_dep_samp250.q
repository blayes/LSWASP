#!/bin/bash
#$ -N dep150_samp250
#$ -q UI-GPU
#$ -pe smp 8
#$ -l h_rt=520:00:00
#$ -l s_rt=520:00:00
#$ -wd /Users/yxu99/dls/code/
#$ -m a
#$ -M yixiang-xu@uiowa.edu
#$ -t 1-1500
#$ -V
#$ -e /Users/yxu99/err/
#$ -o /Users/yxu99/out/
#$ -j y

module load R

R CMD BATCH --no-save --no-restore "--args 4 $SGE_TASK_ID" submit_samp250.R samp/dep_k150_samp250_$SGE_TASK_ID.rout

