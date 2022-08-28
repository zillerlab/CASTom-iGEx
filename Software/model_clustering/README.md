 # Cases stratification based on imputed gene expression
CASTom-iGEx (Module 3) is a pipeline (R based) that uses gene-level T-scores scaled by their association with trait of interest to cluster patients using a graph-based clustering technique. These groups are then tested for association with endophenotype differences and genes/pathway scores as well as group-specific treatment response when data available. If no clinical information are present, plausible differences in endophenotypes are detected with the approximation of gene risk-scores. Two versions are available dependening of data structure (single or multiple cohorts). 

## Input Files
- **Genotype matrix** (*--genoDat_file*): dosages for each chromosome (compressed txt) without variants name/position (variants x samples).  *NOTE: SNPs must match with the train genotype data, file must end with chr<>_matrix.txt.gz*
- **Phenotype matrix**: columns must contain `Individual_ID` plus any phenotype to test the association (phenotypes + 1 x samples). This matrix can include multiple phenotypes to be tested. 
- **Phenotype description**: rows refers to phenotypes to be tested. Columns must include: `pheno_id`, `FieldID`, `Field`,  `transformed_type`;  `pheno_id` is used to match columns name in Phenotype matrix, transformed_type is a charachter defining the type of data (continous, binary ecc.)
- **Covariate matrix** (*--covDat_file*): covariates to correct for in the association analysis (covariats + IDs x samples). Columns must contain `Individual_ID` and `genoSample_ID` to match genotype plus covariates to correct for in the phenotype association. Column `Dx` (0 control 1 case) is optional, if present is used to build the reference set when computing T-scores. *Note: samples in genotype and phenotype matrix are matched based on covariate matrix*
- **Reactome Pathway annotation** (*--reactome_file*): .gmt file can be downloaded from https://reactome.org/download-data/ (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **GO Pathway annotation** (*--GOterms_file*): .RData file, can be obtained using *Annotate_GOterm_run.R*, each pathway is a entry in the list with `GOID` `Term` `Ontology` `geneIds` elemnets (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **Custom pathway**: .RData file, similar to GO structure, each pathway is a list entry with `name` and `geneIds` elements. Available for WikiPathways (2019) in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/)

## Workflow (single cohort)
### Clustering
cluster_PGmethod_corrPCs_run.R

### associate clustering with endophenotype:
cluster_associatePhenoGLM_run.R

### associate clustering with molecular features (genes/pathwayScores):
It initially corrects for PCs, uses wilcoxon test and combined in loci/macrogroups

- cluster_associateFeat_corrPCs_run.R
- filter_pathway_jaccard_sim_run.R (filter pathways based on overlap and min/max number of genes)
- cluster_associatePath_corrPCs_run.R (merge GO and Reactome)

### treatment response
cluster_treatmentResponseAnalysis_run.R

### Predict on external cohort
Corrects for new cohort PCs before projecting clustering

cluster_PGmethod_corrPCs_predict_run.R

Evaluate cluster on external cohort: additional phenotype, percentage of repr, n. of loci

cluster_predict_evaluate_run.R

Endophenotype difference for a predicted cluster, phenoInfo not available

cluster_predict_associatePhenoGLM_run.R


## Gene-RS computation and clustering differences
- compute_risk_score_corrPCs_run.R
- evaluate_risk_score_run.R
- plot_evaluate_risk_score_run.R (Figures/)
- cluster_associatePhenoGLM_run.R/cluster_associatePhenoGLM_multipleCohorts_metaAnalysis_run.R (as before)
- compare_endophenotypeAnalysis_clusterRiskScore_run.R (Figures/)
- plot_precision_risk_score_groupSpec_run.R (Figures/)

***
***

## Workflow (multiple cohorts)
### Detect outliers
- detect_outliers_corrPCs_multipleCohorts_run.R
- combine_outliers_cluster_run.R

### Clustering
cluster_PGmethod_corrPCs_multipleCohorts_run.R

### Predict on external cohort
cluster_PGmethod_corrPCs_predict_run.R

### associate clustering with molecular features (genes/pathwayScores):
- cluster_associateFeat_corrPCs_multipleCohorts_run.R
- filter_pathway_jaccard_sim_run.R (filter pathways based on overlap and min/max number of genes)
- cluster_associatePath_corrPCs_multipleCohort_run.R (merge GO and Reactome)

### compute gene-RS 
- compute_risk_score_corrPCs_multipleCohorts_run.R (across all cohort together)
- cluster_associatePhenoGLM_run.R (evaluate group differences)



