#!/bin/bash
#$ -N sim_accuracy_samp500
#$ -q all.q
#$ -l h_rt=520:00:00
#$ -l s_rt=520:00:00
#$ -wd /Users/yxu99/dls/code/
#$ -m a
#$ -M yixiang-xu@uiowa.edu
#$ -t 1-3
#$ -V
#$ -e /Users/yxu99/err/
#$ -o /Users/yxu99/out/
#$ -j y

module load R

R CMD BATCH --no-save --no-restore "--args 7 $SGE_TASK_ID" submit_samp500.R samp/sim_accur_samp500_$SGE_TASK_ID.rout

