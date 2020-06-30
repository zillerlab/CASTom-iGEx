#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/Tscore_GTEx_targetData_%x_splitGenes%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/Tscore_GTEx_targetData_%x_splitGenes%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=3000MB
#SBATCH --cpus-per-task=30

module load R

cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/

id_t=$1

readarray -t tissues < OUTPUT_GTEx/Tissue_CADgwas
t=$(eval echo "\${tissues[${id_t}-1]}")

echo ${t}

inputfile=OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/
mkdir -p OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_subset/

inputfile_list=()

for i in $(seq 100)
do
	inputfile_list+=(${inputfile}split${i}_predictedExpression_filt.txt)
done

Rscript  ${git_fold}/Tscore_splitGenes_targetData_run.R --input_file ${inputfile_list[@]} --nFolds 10 --perc_comp 0.7 --ncores 30 --covDat_file INPUT_DATA_GTEx/CAD/Covariates/UKBB/covariateMatrix_latestW.txt --covDat_notref_file INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_clustering/covariateMatrix_CADsubset.txt --outFold OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_subset/ --split_gene_id ${SLURM_ARRAY_TASK_ID} --split_tot 100


