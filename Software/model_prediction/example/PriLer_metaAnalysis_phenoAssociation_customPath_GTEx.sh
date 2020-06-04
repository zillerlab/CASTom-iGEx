#!/bin/bash
#SBATCH -o /home/luciat/eQTL_PROJECT/err_out_fold/meta_analysis_customPath_%x_GTEx.out
#SBATCH -e /home/luciat/eQTL_PROJECT/err_out_fold/meta_analysis_customPath_%x_GTEx.err
#SBATCH -N 1
#SBATCH --mem=20G
#SBATCH -t 48:00:00

module load pre2019 2019
module load R

cd /home/luciat/eQTL_PROJECT/

readarray -t cohorts < INPUT_DATA/SCZ_cohort_names
name_cohorts=${cohorts[@]}

git_fold=/home/luciat/priler_project/Software/model_prediction/
file_path_name=$1
path_name=$2
tissue=$3
fold_out=OUTPUT_GTEx/predict_PGC/${tissue}/200kb/PGC_GWAS_bin1e-2/Meta_Analysis_SCZ/devgeno0.01_testdevgeno0/
fold_in=OUTPUT_GTEx/predict_PGC/${tissue}/200kb/PGC_GWAS_bin1e-2/${c}/devgeno0.01_testdevgeno0/

file_res=()
file_pheno=()

for c in ${name_cohorts[@]}
do
	echo ${c}
	file_res+=(${fold_in}pval_SCZ_pheno_covCorr_customPath_${path_name}.RData)
	file_pheno+=(INPUT_DATA/Covariates/${c}.phenoMatrix_old.txt)
done

# correct for covariates
Rscript ${git_fold}pheno_association_customPath_metaAnalysis_run.R --res_cohorts ${file_res[@]} --phenoDatFile_cohorts ${file_pheno[@]} --phenoName Dx --outFold ${fold_out} --pathwayStructure_file /home/luciat/priler_project/refData/${file_path_name} --geneSetName ${path_name}  --cov_corr T --name_cohort ${name_cohorts[@]}


