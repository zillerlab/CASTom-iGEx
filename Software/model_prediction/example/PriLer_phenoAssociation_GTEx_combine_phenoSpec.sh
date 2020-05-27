#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/phenoAssociation_combine_t%a_%x.out
#SBATCH -e /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/err_out_fold/phenoAssociation_combine_t%a_%x.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=1


module load R

cd /psycl/g/mpsziller/lucia/UKBB/eQTL_PROJECT/

id_t=${SLURM_ARRAY_TASK_ID}
readarray -t tissues < OUTPUT_GTEx/Tissue_noGWAS
t=$(eval echo "\${tissues[${id_t}-1]}")
name_pheno=$1

fold=OUTPUT_GTEx/predict_UKBB/${t}/200kb/noGWAS/devgeno0.01_testdevgeno0/${name_pheno}_pheno/
ref_fold=/psycl/g/mpsziller/lucia/refData/
git_fold=/psycl/g/mpsziller/lucia/priler_project/Software/model_prediction/

# correct for covariates
Rscript ${git_fold}pheno_association_combine_largeData_run.R --names_file ${name_pheno} --tscoreFold ${fold}Association_tscore_res/ --pathScoreFold_Reactome ${fold}/Association_reactome_res/ --pathScoreFold_GO ${fold}Association_GO_res/ --outFold ${fold} --cov_corr T --phenoAnn_file INPUT_DATA/Covariates/phenotypeDescription_${name_pheno}.txt --reactome_file ${ref_fold}ReactomePathways.gmt --GOterms_file ${ref_fold}GOterm_geneAnnotation_allOntologies.RData



