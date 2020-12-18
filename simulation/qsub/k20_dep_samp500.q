#!/bin/bash
#$ -N dep20_samp500
#$ -q UI-GPU
#$ -pe smp 8
#$ -l h_rt=520:00:00
#$ -l s_rt=520:00:00
#$ -wd /Users/yxu99/dls/code/
#$ -m a
#$ -M yixiang-xu@uiowa.edu
#$ -t 1-200
#$ -V
#$ -e /Users/yxu99/err/
#$ -o /Users/yxu99/out/
#$ -j y

module load R

R CMD BATCH --no-save --no-restore "--args 1 $SGE_TASK_ID" submit_samp500.R samp/dep_k20_samp500_$SGE_TASK_ID.rout

