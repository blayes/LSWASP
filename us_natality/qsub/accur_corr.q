#!/bin/bash
#$ -N real_accuracy_corr
#$ -q all.q
#$ -pe smp 8
#$ -l h_rt=520:00:00
#$ -l s_rt=520:00:00
#$ -wd /Users/yxu99/real_dls_overlap/code/
#$ -m a
#$ -M yixiang-xu@uiowa.edu
#$ -t 1-2
#$ -V
#$ -e /Users/yxu99/err/
#$ -o /Users/yxu99/out/
#$ -j y

module load R

R CMD BATCH --no-save --no-restore "--args 15 $SGE_TASK_ID" submit.R samp/real_accur_corr_$SGE_TASK_ID.rout

