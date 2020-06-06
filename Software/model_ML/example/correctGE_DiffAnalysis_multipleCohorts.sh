#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/corrGE_diffAnalysis_target1235_t%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/corrGE_diffAnalysis_target1235_t%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=1

module load R

cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/

git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_ML/
ref_fold=/psycl/g/mpsziller/lucia/priler_project/refData/

id_t=${SLURM_ARRAY_TASK_ID}
readarray -t tissues < OUTPUT_GTEx/Tissue_CADgwas
t=$(eval echo "\${tissues[${id_t}-1]}")

fold_cov=INPUT_DATA_GTEx/CAD/Covariates/
fold=OUTPUT_GTEx/predict_CAD/${t}/200kb/CAD_GWAS_bin5e-2/Meta_Analysis_CAD/devgeno0.01_testdevgeno0/

name_cohorts=(CG German1 German2 German3 German4 German5 LURIC MG WTCCC)
sample_ref_f=()
for c in ${name_cohorts[@]}
do
	sample_ref_f+=(${fold_cov}${c}/covariateMatrix_base_seed1235.txt)
done

Rscript ${git_fold}PathwayDiffAnalysis_correctedGeneExp_multipleCohorts_run.R --name_cohorts ${name_cohorts[@]} --sampleAnn_file ${fold_cov}Meta_Analysis_CAD/covariateMatrix.txt --sampleAnn_ref_file ${sample_ref_f[@]} --input_file ${fold}corrected_predictedExpression.txt.gz --outFold ${fold}target1235_ --reactome_file ${ref_fold}ReactomePathways.gmt --GOterms_file ${ref_fold}GOterm_geneAnnotation_allOntologies.RData 
