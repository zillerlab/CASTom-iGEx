#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/pheno_pathGO_%x_split%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/pheno_pathGO_%x_split%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2

module load R
cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/
git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_prediction/

name_pheno=$1
mkdir -p OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/${name_pheno}_pheno/Association_GO_res/

fold=OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/${name_pheno}_pheno/

# correct for covariates
Rscript ${git_fold}pheno_association_pathscore_largeData_run.R --inputInfoFile ${fold}pathScore_GO_info.RData --covDat_file INPUT_DATA/Covariates/covariateMatrix_${name_pheno}.txt --phenoDat_file INPUT_DATA/Covariates/phenotypeMatrix_${name_pheno}.txt --names_file ${name_pheno} --inputFile ${fold}Pathway_GO_scores_splitPath${SLURM_ARRAY_TASK_ID}.RData --outFile ${fold}Association_GO_res/pval_pathScore_GO_splitPath${SLURM_ARRAY_TASK_ID}_ --cov_corr T --sampleAnn_file INPUT_DATA/Covariates/covariateMatrix_${name_pheno}.txt --ncores 2 --phenoAnn_file INPUT_DATA/Covariates/phenotypeDescription_${name_pheno}.txt  --functR ${git_fold}pheno_association_functions.R --path_type GO  

