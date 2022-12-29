 # Cases stratification based on imputed gene expression
CASTom-iGEx (Module 3) is a pipeline (R based) that uses gene-level T-scores, corrected for PCs and scaled by their association with trait of interest, to cluster patients using a graph-based clustering technique. These groups are then tested for association with endophenotype differences and genes/pathway scores as well as group-specific treatment response when data available. If no clinical information are present, plausible differences in endophenotypes are detected with the approximation of gene risk-scores. Two versions are available dependening of data structure (single or multiple cohorts). 

## Input 
### Data
- **Sample matrix** (*--sampleAnnFile*): .txt file tab separated, includes the subset of samples to be clustered. Columns must contain `Individual_ID` and `Dx` that refers to Cases (`Dx=1`) and Controls (`Dx=0`).
- **Input matrix** (*--inputFile*): `predictedTscores.txt` tab separated (small n. samples) OR `predictedTscore_splitGenes` common name of split predicted gene T-scores in .RData format (large n. samples). These files are produced from `Tscore_PathScore_diff_run.R` or `Tscore_splitGenes_run.R` in [CASTom-iGEx Module 2](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_prediction).  
- **TWAS and PALAS results** (*--pvalresFile*):

### Parameters
- **type of clustering** (*--type_cluster*): Indicates if all samples in *--sampleAnnFile* are clustered (`All`), only affected individuals (`Cases`) or only non-affected individuals (`Controls`).
- **type of similarity** (*--type_sim*): if `HK`, heat kernel similarity is used. If `ED`, opposite of euclidian distance is used as similarity.
- **type of data** (*--type_data*):
- **type of input** (*--type_input*): how input should be processed
- **K nearest neighbour** (*--kNN_par*): number of nearest neighbours considered to compute both heat kernel similarity and shared nearest neighbour similarity.


FROM HERE

- **Phenotype matrix**: columns must contain `Individual_ID` plus any phenotype to test the association (phenotypes + 1 x samples). This matrix can include multiple phenotypes to be tested. 
- **Phenotype description**: rows refers to phenotypes to be tested. Columns must include: `pheno_id`, `FieldID`, `Field`,  `transformed_type`;  `pheno_id` is used to match columns name in Phenotype matrix, transformed_type is a charachter defining the type of data (continous, binary ecc.)
- **Covariate matrix** (*--covDat_file*): covariates to correct for in the association analysis (covariats + IDs x samples). Columns must contain `Individual_ID` and `genoSample_ID` to match genotype plus covariates to correct for in the phenotype association. Column `Dx` (0 control 1 case) is optional, if present is used to build the reference set when computing T-scores. *Note: samples in genotype and phenotype matrix are matched based on covariate matrix*
- **Reactome Pathway annotation** (*--reactome_file*): .gmt file can be downloaded from https://reactome.org/download-data/ (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **GO Pathway annotation** (*--GOterms_file*): .RData file, can be obtained using *Annotate_GOterm_run.R*, each pathway is a entry in the list with `GOID` `Term` `Ontology` `geneIds` elemnets (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **Custom pathway**: .RData file, similar to GO structure, each pathway is a list entry with `name` and `geneIds` elements. Available for WikiPathways (2019) in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/)

### Optional: Initial filtering if datasets are not harmonized
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
- tissue_name_filter_genes_matched_datasets.txt 

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
- tissue_name_filter_path_type_path_matched_datasets.txt

### Optional: Clustering based on genetic principal componenets
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


## Single cohort (tissue-specific)
### 1) Clustering
Cluster individuals based on tissue-specific imputed gene expression. This script intially clump genes based on correlation, standardize each genes, correct for PCs and rescaled them by TWAS disease-specifc Z-statistic.
- *--covDatFile*
- *--split_tot*
- *--pval_id*
- *--corr_thr*
- *--min_genes_path*
- *--exclude_MHC*
- *--capped_zscore*
- *--geneRegionFile*

```sh
./cluster_PGmethod_corrPCs_run.R \
    --inputFile \
    --sampleAnnFile \
    --tissues_name \
    --covDatFile \
    --type_cluster \
    --split_tot (default 0) \
    --pvalresFile \
    --pval_id \
    --pval_thr (default 1) \
    --corr_thr (default -1) \
    --functR ./clustering_functions.R \
    --type_data help = "tscore, path_Reactome or path_GO" \
    --type_sim \
    --type_input "original or zscaled" \
    --kNN_par \
    --min_genes_path (default 1) \
    --exclude_MHC (default FALSE) "if true, MHC region excluded (only ossible for tscore)"\
    --capped_zscore (default FALSE) "if true, zstat is capped based on distribution"\
    --geneRegionFile (default NULL) "used if tscore and exclude_MHC"\
    --outFold
```

### Optional: Clustering without PCs correction
Cluster as before but without correcting for PCs, this script is used to benchmark results and compare them with the clustering correcting for PCs. Inputs specifics as previous script.

```sh
./cluster_PGmethod_run.R \
    --inputFile \
    --sampleAnnFile \
    --tissues_name \
    --covDatFile \
    --type_cluster \
    --split_tot (default 0) \
    --pvalresFile \
    --pval_id \
    --pval_thr (default 1) \
    --corr_thr (default -1) \
    --functR ./clustering_functions.R \
    --type_data \
    --type_sim \
    --type_input \
    --kNN_par \
    --min_genes_path (default 1) \
    --exclude_MHC (default FALSE) \
    --capped_zscore (default FALSE) \
    --geneRegionFile (default NULL) \
    --outFold
```


### 2.1) Associate clusters with molecular features (genes/pathwayScores):
It initially corrects for PCs, uses wilcoxon test and combined in loci/macrogroups

- cluster_associateFeat_corrPCs_run.R
- filter_pathway_jaccard_sim_run.R (filter pathways based on overlap and min/max number of genes)
- cluster_associatePath_corrPCs_run.R (merge GO and Reactome)

### 2.2) Associate clusters with endophenotypes
cluster_associatePhenoGLM_run.R

### 2.3) Find differential treatment response among clusters
cluster_treatmentResponseAnalysis_run.R

### 2.4) Drug repositiong based on cluster-specific pathways
pathSEA_path_group_run.R

### 3) Project on external cohorts
Corrects for new cohort PCs before projecting clustering

cluster_PGmethod_corrPCs_predict_run.R

Evaluate cluster on external cohort: additional phenotype, percentage of repr, n. of loci

cluster_predict_evaluate_run.R

Endophenotype difference for a predicted cluster, phenoInfo not available

cluster_predict_associatePhenoGLM_run.R

### 4) Compute gene-risk score
compute_risk_score_corrPCs_run.R
***

## Optional: Evaluate gene-risk score 
- evaluate_risk_score_run.R
- plot_evaluate_risk_score_run.R (Figures/)
- cluster_associatePhenoGLM_run.R/cluster_associatePhenoGLM_multipleCohorts_metaAnalysis_run.R (as before)
- compare_endophenotypeAnalysis_clusterRiskScore_run.R (Figures/)
- plot_precision_risk_score_groupSpec_run.R (Figures/)

***

## Multiple cohorts (tissue-specific)
### 0) Detect outliers
- detect_outliers_corrPCs_multipleCohorts_run.R
- combine_outliers_cluster_run.R

### 1) Clustering
cluster_PGmethod_corrPCs_multipleCohorts_run.R

### 2.1) Associate clusters with molecular features (genes/pathwayScores):
- cluster_associateFeat_corrPCs_multipleCohorts_run.R
- filter_pathway_jaccard_sim_run.R (filter pathways based on overlap and min/max number of genes)
- cluster_associatePath_corrPCs_multipleCohort_run.R (merge GO and Reactome)

### 2.2) Associate clusters with endophenotypes
cluster_associatePhenoGLM_run.R

### 3) Project on external cohorts
cluster_PGmethod_corrPCs_predict_run.R

### 4) Compute gene-risk score
- compute_risk_score_corrPCs_multipleCohorts_run.R (across all cohort together)
- cluster_associatePhenoGLM_run.R (evaluate group differences)



