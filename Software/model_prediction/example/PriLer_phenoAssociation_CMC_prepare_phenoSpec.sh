#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/phenoAssociation_prepare_CMC_%x.out
#SBATCH -e /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/phenoAssociation_prepare_CMC_%x.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=60G
#SBATCH --cpus-per-task=1

module load R
cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/

name_pheno=$1

git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_prediction/ 
fold_train=OUTPUT_CMC/train_CMC/200kb/
fold=OUTPUT_CMC/predict_UKBB/200kb/devgeno0.01_testdevgeno0/${name_pheno}_pheno/
ref_fold=/psycl/g/mpsziller/lucia/refData/

Rscript ${git_fold}pheno_association_prepare_largeData_run.R --split_tot 100 --geneAnn_file ${fold_train}resPrior_regEval_allchr.txt --inputFold ${fold} --outFold ${fold} --GOterms_file ${ref_fold}GOterm_geneAnnotation_allOntologies.RData --reactome_file ${ref_fold}ReactomePathways.gmt  --sampleAnn_file INPUT_DATA/Covariates/covariateMatrix_${name_pheno}.txt 
