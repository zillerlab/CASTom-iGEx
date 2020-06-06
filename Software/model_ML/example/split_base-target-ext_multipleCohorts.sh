#!/bin/bash
#SBATCH -o /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/split_base_target_ext_c%a.out
#SBATCH -e /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/err_out_fold/split_base_target_ext_c%a.err
#SBATCH --time=7-0
#SBATCH --nodes=1
#SBATCH --mem=10G

module load R
cd /psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/

readarray -t cohorts < INPUT_DATA_GTEx/CAD/CAD_cohort_name
id_c=${SLURM_ARRAY_TASK_ID}
cohort=$(eval echo "\${cohorts[${id_c}-1]}")

fold_git=/psycl/g/mpsziller/lucia/priler_project/Software/model_ML/

Rscript ${fold_git}split_base-target-ext_multipleCohorts_run.R --covDat_file INPUT_DATA_GTEx/CAD/Covariates/${cohort}/covariateMatrix.txt --perc_base 0.7 --perc_ext 0.2 --seed_fixed 1235 --outFold INPUT_DATA_GTEx/CAD/Covariates/${cohort}/

