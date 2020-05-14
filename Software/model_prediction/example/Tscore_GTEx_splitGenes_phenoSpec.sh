#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/Tscore_GTEx_%x_splitGenes%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/Tscore_GTEx_%x_splitGenes%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=3000MB
#SBATCH --cpus-per-task=30

module load R

cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/

id_t=$1
name_pheno=$2
perc_comp=$3

readarray -t tissues < OUTPUT_GTEx/Tissue_noGWAS
t=$(eval echo "\${tissues[${id_t}-1]}")

echo ${t}

mkdir -p OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/${name_pheno}_pheno/

inputfile=OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/


inputfile_list=()

for i in $(seq 100)
do
	inputfile_list+=(${inputfile}split${i}_predictedExpression_filt.txt)
done

Rscript  RSCRIPTS/SCRIPTS_v2/Tscore_splitGenes_run.R --input_file ${inputfile_list[@]} --nFolds 10 --perc_comp ${perc_comp} --ncores 30 --covDat_file INPUT_DATA/Covariates/covariateMatrix_${name_pheno}.txt --outFold OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/${name_pheno}_pheno/ --split_gene_id ${SLURM_ARRAY_TASK_ID} --split_tot 100


