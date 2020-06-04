#!/bin/bash
#SBATCH -o /home/luciat/eQTL_PROJECT/err_out_fold/phenoAssociation_customPath_%x_c%a_GTEx.out
#SBATCH -e /home/luciat/eQTL_PROJECT/err_out_fold/phenoAssociation_customPath_%x_c%a_GTEx.err
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH -t 24:00:00

module load pre2019 2019
module load R

cd /home/luciat/eQTL_PROJECT/

id_c=${SLURM_ARRAY_TASK_ID}
readarray -t cohorts < INPUT_DATA/SCZ_cohort_names
c=$(eval echo "\${cohorts[${id_c}-1]}")

git_fold=/home/luciat/priler_project/Software/model_prediction/
file_path_name=$1
path_name=$2
tissue=$3
fold=OUTPUT_GTEx/predict_PGC/${tissue}/200kb/PGC_GWAS_bin1e-2/${c}/devgeno0.01_testdevgeno0/
fold_train=OUTPUT_GTEx/train_GTEx/${tissue}/200kb/PGC_GWAS_bin1e-2/

Rscript ${git_fold}pheno_association_smallData_customPath_run.R --covDat_file INPUT_DATA/Covariates/${c}.covariateMatrix_old.txt --phenoDat_file INPUT_DATA/Covariates/${c}.phenoMatrix_old.txt --geneAnn_file ${fold_train}resPrior_regEval_allchr.txt --inputFold ${fold} --outFold ${fold} --pathwayStructure_file /home/luciat/priler_project/refData/${file_path_name} --cov_corr T --functR ${git_fold}pheno_association_functions.R --sampleAnn_file INPUT_DATA/Covariates/${c}.covariateMatrix_old.txt --names_file SCZ_pheno --phenoAnn_file INPUT_DATA/Covariates/phenotypeDescription_PGCcohorts.csv --geneSetName ${path_name} 


