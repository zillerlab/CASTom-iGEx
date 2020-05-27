#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/split_data_base_target_%x.out
#SBATCH -e /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/split_data_base_target_%x.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=10G

module load R
cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/

name_pheno=$1
perc_comp=$2
t=$3

cov_fold=INPUT_DATA/Covariates/
git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_ML/
input_file=()
for i in $(seq 100)
do
input_file+=(OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/split${i}_predictedExpression_filt.txt)
done

Rscript ${git_fold}split_base-target-ext_run.R --input_file ${input_file[@]} --covDat_file ${cov_fold}covariateMatrix_${name_pheno}.txt --nFolds 10 --perc_comp ${perc_comp} --seed_fixed 1234 --thr_n 10000 --outFold ${cov_fold}${name_pheno}_

