#!/bin/bash
#$ -N ind50
#$ -q UI-GPU
#$ -pe smp 8
#$ -l h_rt=520:00:00
#$ -l s_rt=520:00:00
#$ -wd /Users/yxu99/dls/code/
#$ -m a
#$ -M yixiang-xu@uiowa.edu
#$ -t 1-500
#$ -V
#$ -e /Users/yxu99/err/
#$ -o /Users/yxu99/out/
#$ -j y

module load R

R CMD BATCH --no-save --no-restore "--args 3 $SGE_TASK_ID" submit.R samp/ind_k50_$SGE_TASK_ID.rout

