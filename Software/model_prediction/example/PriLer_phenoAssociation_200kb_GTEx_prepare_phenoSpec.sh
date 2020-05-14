#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/phenoAssociation_prepare_t%a_%x.out
#SBATCH -e /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/phenoAssociation_prepare_t%a_%x.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=60G
#SBATCH --cpus-per-task=1

module load R
cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/


readarray -t tissues < OUTPUT_GTEx/Tissue_noGWAS
t=$(eval echo "\${tissues[${id_t}-1]}")
name_pheno=$1
 
fold_train=OUTPUT_GTEx/train_GTEx/${t}/200kb/noGWAS/
fold=OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/${name_pheno}_pheno/
ref_fold=/psycl/g/mpsziller/lucia/refData/

Rscript RSCRIPTS/SCRIPTS_v2/pheno_association_prepare_run.R --split_tot 100 --geneAnn_file ${fold_train}resPrior_regEval_allchr.txt --inputFold ${fold} --outFold ${name_pheno} --GOterms_file ${ref_fold}GOterm_geneAnnotation_allOntologies.RData --reactome_file ${ref_fold}ReactomePathways.gmt  --sampleAnn_file INPUT_DATA/Covariates/covariateMatrix_${name_pheno}.txt 
