#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/PathScore_GTEx_%x.out
#SBATCH -e /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/PathScore_GTEx_%x.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=35G
#SBATCH --cpus-per-task=10

module load R

id_t=$1
name_pheno=$2

cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/

readarray -t tissues < OUTPUT_GTEx/Tissue_noGWAS
t=$(eval echo "\${tissues[${id_t}-1]}")

fold_ann=/psycl/g/mpsziller/lucia/refData/
inputfold=OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/${name_pheno}_pheno/
covfold=INPUT_DATA/Covariates/


Rscript  RSCRIPTS/SCRIPTS_v2/PathwayScores_splitGenes_run.R --ncores 10 --input_file ${inputfold}predictedTscores_splitGenes  --covDat_file ${covfold}covariateMatrix_${name_pheno}.txt  --outFold ${inputfold}  --split_tot 100 --reactome_file ${fold_ann}ReactomePathways.gmt --GOterms_file ${fold_ann}GOterm_geneAnnotation_allOntologies.RData


