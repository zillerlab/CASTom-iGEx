#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/corr_predExpr_CADcohorts_t%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/corr_predExpr_CADcohorts_t%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1

module load R

cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_prediction/

id_t=${SLURM_ARRAY_TASK_ID}
readarray -t tissues < OUTPUT_GTEx/Tissue_CADgwas
t=$(eval echo "\${tissues[${id_t}-1]}")

fold_cov=INPUT_DATA_GTEx/CAD/Covariates/
fold_input=OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/

name_cohorts=(CG German1 German2 German3 German4 German5 LURIC MG WTCCC)
input_f=()
cov_f=()
sample_f=()
for c in ${name_cohorts[@]}
do
	sample_f+=(${fold_cov}${c}/covariateMatrix.txt)
	cov_f+=(${fold_cov}${c}/covariateMatrix.txt)
	input_f+=(${fold_input}${c}/predictedExpression.txt.gz)
done

fold_output=${fold_input}Meta_Analysis_CAD/devgeno0.01_testdevgeno0/
fold_output_ann=${fold_cov}Meta_Analysis_CAD/
mkdir -p ${fold_output}


Rscript ${git_fold}correct_predictedGeneExpression_multipleCohorts_run.R --name_cohorts ${name_cohorts[@]} --seed_umap 50 --n_comp_umap 2 --covDat_file ${cov_f[@]} --sampleAnn_file ${sample_f[@]} --input_file ${input_f[@]} --outFold ${fold_output} --outFoldAnn ${fold_output_ann}
