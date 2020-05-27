#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/pheno_tscore_%x_split%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/pheno_tscore_%x_split%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2

module load R
cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/
git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_prediction/

name_pheno=$1
id_t=$2

readarray -t tissues < OUTPUT_GTEx/Tissue_noGWAS
t=$(eval echo "\${tissues[${id_t}-1]}")
mkdir -p OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/${name_pheno}_pheno/Association_tscore_res/

fold=OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/${name_pheno}_pheno/

# correct for covariates
Rscript ${git_fold}pheno_association_tscore_largeData_run.R --split_tot 100 --inputInfoFile ${fold}tscore_info.RData --covDat_file INPUT_DATA/Covariates/covariateMatrix_${name_pheno}.txt --phenoDat_file INPUT_DATA/Covariates/phenotypeMatrix_${name_pheno}.txt --names_file ${name_pheno}  --inputFile ${fold}predictedTscores_splitGenes${SLURM_ARRAY_TASK_ID}.RData --outFile ${fold}Association_tscore_res/pval_tscore_splitGene${SLURM_ARRAY_TASK_ID}_ --split_gene_id ${SLURM_ARRAY_TASK_ID} --cov_corr T --sampleAnn_file INPUT_DATA/Covariates/covariateMatrix_${name_pheno}.txt --ncores 2 --phenoAnn_file INPUT_DATA/Covariates/phenotypeDescription_${name_pheno}.txt  --functR ${git_fold}pheno_association_functions.R  

