 # Cases stratification based on imputed gene expression
CASTom-iGEx (Module 3) is a pipeline (R based) that uses gene-level T-scores, corrected for PCs and scaled by their association with trait of interest, to cluster patients using a graph-based clustering technique. These groups are then tested for association with endophenotype differences and genes/pathway scores as well as group-specific treatment response when data available. If no clinical information are present, plausible differences in endophenotypes are detected with the approximation of gene risk-scores. Two versions are available dependening of data structure (single or multiple cohorts). 

## Input Files
- **Sample matrix** (*--sampleAnnFile*): .txt file tab separated, includes the subset of samples to be clustered. Columns must contain `Individual_ID` and `Dx` that refers to Cases (`Dx=1`) and Controls (`Dx=0`).

## Input parameters
- **type of clustering** (*--type_cluster*): Indicates if all samples in *--sampleAnnFile* are clustered (`All`), only affected individuals (`Cases`) or only non-affected individuals (`Controls`).
- **type of similarity** (*--type_sim*): if `HK`, heat kernel similarity is used. If `ED`, opposite of euclidian distance is used as similarity.
- **K nearest neighbour** (*--kNN_par*): number of nearest neighbours considered to compute both heat kernel similarity and shared nearest neighbour similarity.


FROM HERE

- **Phenotype matrix**: columns must contain `Individual_ID` plus any phenotype to test the association (phenotypes + 1 x samples). This matrix can include multiple phenotypes to be tested. 
- **Phenotype description**: rows refers to phenotypes to be tested. Columns must include: `pheno_id`, `FieldID`, `Field`,  `transformed_type`;  `pheno_id` is used to match columns name in Phenotype matrix, transformed_type is a charachter defining the type of data (continous, binary ecc.)
- **Covariate matrix** (*--covDat_file*): covariates to correct for in the association analysis (covariats + IDs x samples). Columns must contain `Individual_ID` and `genoSample_ID` to match genotype plus covariates to correct for in the phenotype association. Column `Dx` (0 control 1 case) is optional, if present is used to build the reference set when computing T-scores. *Note: samples in genotype and phenotype matrix are matched based on covariate matrix*
- **Reactome Pathway annotation** (*--reactome_file*): .gmt file can be downloaded from https://reactome.org/download-data/ (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **GO Pathway annotation** (*--GOterms_file*): .RData file, can be obtained using *Annotate_GOterm_run.R*, each pathway is a entry in the list with `GOID` `Term` `Ontology` `geneIds` elemnets (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **Custom pathway**: .RData file, similar to GO structure, each pathway is a list entry with `name` and `geneIds` elements. Available for WikiPathways (2019) in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/)

### Initial filtering if datasets are not harmonized (optional)
The following two scripts are used when gene risk scores are predicted on a data set (e.g. PGC) but TWAS estimates are obtained from another data set (e.g. UKBB) and the two initial PriLer models were not harmonized per SNPs. This step is prior the clustering computation that should be based on a subset of correlated genes.

#### Compare imputed genes
Compute genes correlation imputed from 2 different models. Genes are imputed on the reference panel from which the gene expression models are estimated, see [CASTom-iGEx Module 1](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_training)
- *--geneExpPred_file*: predicted gene expression on the same reference panel (e.g. GTEx) for a certain tissue but based on 2 different harmonization (reference panel + data set 1 and reference panel + data set 2). Complete file path for the two imputations.

```sh
./compare_geneExp_matchedDataset_run.R \
    --geneExpPred_file (2 files necessary) \
    --corr_thr (default 0.8) \
    --tissue_name \
    --outFold 
```
The output includes (saved in *--outFold*):
- *tissue_name*_filter_genes_matched_datasets.txt 

#### Compare imputed pathways
Compute pathways correlation imputed from 2 different models. Pathways-scores are computed on the reference panel from which the gene expression models are estimated (see [CASTom-iGEx Module 1](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_training)), after the computation of gene T-scores (see [CASTom-iGEx Module 2](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_prediction))
- *--pathScore_file*: predicted pathway-score on the same reference panel (e.g. GTEx) for a certain tissue but based on 2 different harmonization (reference panel + data set 1 and reference panel + data set 2). Complete file path for the two imputations.
- *--type_path*: pathway database name, Reactome or GO.

```sh
./compare_pathScore_matchedDataset_run.R \
    --pathScore_file (2 files necessary) \
    --corr_thr (default 0.8) \
    --type_path \
    --tissue_name \
    --outFold 
```
The output includes (saved in *--outFold*):
- *tissue_name*_filter_path_*type_path*_matched_datasets.txt

### Clustering based on genetic principal componenets (optional)
Cluster individuals based on genetic principal componenets. This script is used to benchmark results and observe the overlap with tissue-specific clustering.
- *--PCs_input_file*: .RData object containing matrix (nsamples x PCs) rownames must be `Individual_ID`.
- *--sampleOutFile*: .txt tab separated file, list of samples to be removed from clustering, columns must include `Individual_ID`.

```sh
./cluster_PGmethod_PCs_run.R \
    --PCs_input_file \
    --sampleAnnFile \
    --sampleOutFile \
    --type_cluster (default "All") \
    --functR ./clustering_functions.R \
    --type_sim (default "HK") \
    --kNN_par (default 30) \
    --outFold
```
The output includes (saved in *--outFold*):
- PCs_cluster*type_cluster*_PGmethod_*type_sim*metric.RData object with the following structure:
    - best_k: single kNN parameter given OR if multiple kNN provided, returns kNN that maximizes Davies-Bouldin index.
    - cl_res: total clustering ouptput including computed similarity matrix.
    - test_cov: chisq-test or kruskal-wallis test of clustering structure and covariates (e.g. PCs or sex).
    - info_tune: metric of clustering given kNN, includes modularity maximized by Louvain clustering
    - feat: names of feautures used for clustering
    - cl_best: optimal (if multiple kNN given) or unique clustering structure obtained
    - Dx_perc: If `--type_cluster All`, calculates the percentages of cases and controls in the found partition
    - samples_id: `Individual_ID` of the considered samples
    - test_diff_gr: kruskal wallis test of each features used to compute the clusters and the clustering structure.
    - gr_input: mean, sd and coefficient of variation (cv) of each feature across clusters.
    - input_data: standardized input features used for clustering


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

### drug repositioning 
pathSEA_path_group_run.R

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



