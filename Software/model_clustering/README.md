 # Cases stratification based on imputed gene expression
CASTom-iGEx (Module 3) is a pipeline (R based) that uses gene-level T-scores, corrected for PCs and scaled by their association with trait of interest, to cluster patients using a graph-based clustering technique. These groups are then tested for association with endophenotype differences and genes/pathway scores as well as group-specific treatment response when data available. If no clinical information are present, plausible differences in endophenotypes are detected with the approximation of gene risk-scores. Two versions are available dependening of data structure (single or multiple cohorts). 

## Input 
### Data
- **Sample matrix** (*--sampleAnnFile*): .txt file tab separated, includes the subset of samples to be clustered. Columns must contain `Individual_ID` and `Dx` that refers to Cases (`Dx=1`) and Controls (`Dx=0`) and principal components. 
- **Input matrix** (*--inputFile*): `predictedTscores.txt` tab separated (small n. samples) OR `predictedTscore_splitGenes` common name of split predicted gene T-scores in .RData format (large n. samples). These files are produced from `Tscore_PathScore_diff_run.R` or `Tscore_splitGenes_run.R` in [CASTom-iGEx Module 2](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_prediction).  
- **TWAS and PALAS results** (*--pvalresFile*): .RData file obtained from `pheno_association_` scripts in [CASTom-iGEx Module 2](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_prediction). Z-statistic from genes or pathways are used to annotate clustering features and re-scale them.
- **Phenotype matrix** (*--phenoDatFile*): columns must contain `Individual_ID` plus any phenotype to test for clustering differences (phenotypes + 1 x samples). This matrix can include multiple phenotypes to be tested. 
- **Phenotype description** (*--phenoDescFile*): rows refers to phenotypes to be tested. Columns must include: `pheno_id`, `FieldID`, `Field`,  `transformed_type`;  `pheno_id` is used to match columns name in Phenotype matrix, transformed_type is a charachter defining the type of data (continous, binary ecc.)

 
### Parameters and auxiliary files
- **type of clustering** (*--type_cluster*): Indicates if all samples in *--sampleAnnFile* are clustered (`All`), only affected individuals (`Cases`) or only non-affected individuals (`Controls`).
- **type of similarity** (*--type_sim*): if `HK`, heat kernel similarity is used. If `ED`, opposite of euclidian distance is used as similarity.
- **type of data** (*--type_data*): indicates the type of input data (*--inputFile*) used to compute the clustering: gene T-score `tscore`, Pathway-score Reactome `path_Reactome` or pathway-score Gene Ontology `path_GO`.
- **type of input** (*--type_input*): indicates how the input should be processed, if `original` the data is not multiplied by Z-statistic nut only standardized, if `zscaled` each feature is first standardized and then multiplied by the correspondig Z-statistic provided in *--pvalresFile*.
- **K nearest neighbour** (*--kNN_par*): number of nearest neighbours considered to compute both heat kernel similarity and shared nearest neighbour similarity.
- **Exclusion of MHC locus** (*--exclude_MHC*): if TRUE and `type_data="tscore"`, exclude genes in MHC locus.
- **Gene coordinate info** (*--geneRegionFile*): used if `type_data="tscore"` and `exclude_MHC=TRUE`. Tab separated file that contains gene location (e.g. resPrior_regEval_allchr.txt from PriLer output [CASTom-iGEx Module 1](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/Software/model_training)).
- **Split of Input Matrix** (*--split_tot*): integer indicating the number of groups the *--inputFile* has been split. Used in a large sample size setting e.g. UK Biobank.
- **Method for community detection** (*--cluster_method*): leiden or louvain, refers to the community detection strategy

***
### Optional: Clustering based on genetic principal componenets
Cluster individuals based on genetic principal componenets. This script is used to benchmark results and observe the overlap with tissue-specific clustering.
- *--PCs_input_file*: .RData object containing matrix (nsamples x PCs) rownames must be `Individual_ID`.
- *--sampleOutFile*: .txt tab separated file, list of samples to be removed from clustering, columns must include `Individual_ID`.

```sh
./cluster_PGmethod_PCs_run.R \
    --PCs_input_file \
    --sampleAnnFile \
    --sampleOutFile \
    --type_cluster \
    --functR ./clustering_functions.R \
    --type_sim (default "HK") \
    --kNN_par (default 20) \
    --cluster_method (default "leiden") \
    --outFold
```
The output includes (saved in *--outFold*):
- PCs_cluster**type_cluster**\_PGmethod\_**type_sim**metric.RData object with the following structure:
    - best_k: single kNN parameter given OR if multiple kNN provided, returns kNN that maximizes Davies-Bouldin index.
    - cl_res: total clustering ouptput including computed similarity matrix.
    - test_cov: chisq-test or kruskal-wallis test of clustering structure and covariates (e.g. PCs or sex).
    - info_tune: metric of clustering given kNN, includes modularity maximized by Louvain/Leiden clustering
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
- *--covDatFile*: optional, additional covariates to test for. NOTE: principal components shoudl be included in *--sampleAnnFile*.
- *--pval_id*: integer indicating the index of the phenotype to be considered in *--pvalresFile* that can include more than one phenotype.
- *--corr_thr*: correlation threshold to clump features. If 1 include all the features.
- *--min_genes_path*: minimum number of genes forming a pathway. Used to filter pathways for clustering with *--type_data* is path_Reactome or path_GO.
- *--capped_zscore*: if TRUE, Z-statistic to rescale features are capped at 0.05% to attenuate extreme associations

```sh
./cluster_PGmethod_corrPCs_run.R \
    --inputFile \
    --sampleAnnFile \
    --tissues_name \
    --covDatFile (default NULL) \
    --type_cluster \
    --split_tot (default 0) \
    --pvalresFile \
    --pval_id (default 1) \
    --pval_thr (default 1) \
    --corr_thr (default -1) \
    --functR ./clustering_functions.R \
    --type_data (default tscore) \
    --type_sim (default HK) \
    --type_input (default original) \
    --kNN_par (default 20) \
    --min_genes_path (default 1) \
    --exclude_MHC (default FALSE) \
    --capped_zscore (default FALSE) \
    --geneRegionFile (default NULL) \
    --cluster_method (default "leiden") \
    --outFold
```
The output includes (saved in *--outFold*):
- **type_data**\_corrPCs\_**type_input**\_cluster**type_cluster**\_PGmethod\_**type_sim**metric.RData object with the following structure:
    - best_k: single kNN parameter given OR if multiple kNN provided, returns kNN that maximizes Davies-Bouldin index.
    - cl_res: total clustering ouptput including computed similarity matrix.
    - test_cov: chisq-test or kruskal-wallis test of clustering structure and covariates (e.g. PCs or sex).
    - info_tune: metric of clustering given kNN, includes modularity maximized by Louvain/Leiden clustering
    - feat: names of feautures used for clustering
    - cl_best: optimal (if multiple kNN given) or unique clustering structure obtained
    - Dx_perc: If `--type_cluster All`, calculates the percentages of cases and controls in the found partition
    - samples_id: `Individual_ID` of the considered samples
    - test_diff_gr: kruskal wallis test of each features used to compute the clusters and the clustering structure.
    - gr_input: mean, sd and coefficient of variation (cv) of each feature across clusters.
    - input_data: standardized input features used for clustering

### Optional: Clustering without PCs correction
Cluster as before but without correcting for PCs, this script is used to benchmark results and compare them with the clustering correcting for PCs. Inputs specifics as previous script.

```sh
./cluster_PGmethod_run.R \
    --inputFile \
    --sampleAnnFile \
    --tissues_name \
    --covDatFile (default NULL) \
    --type_cluster \
    --split_tot (default 0) \
    --pvalresFile \
    --pval_id (default 1) \
    --pval_thr (default 1) \
    --corr_thr (default -1) \
    --functR ./clustering_functions.R \
    --type_data (default "tscore") \
    --type_sim (default "HK") \
    --type_input (default "original") \
    --kNN_par (default 20) \
    --min_genes_path (default 1) \
    --exclude_MHC (default FALSE) \
    --capped_zscore (default FALSE) \
    --geneRegionFile (default NULL) \
    --cluster_method (default "leiden") \
    --outFold
```
The output includes (saved in *--outFold*):
- **type_data**\_**type_input**\_cluster**type_cluster**\_PGmethod\_**type_sim**metric.RData object with the same structure as the previous command.

### 2.1) Associate clusters with molecular features (genes/pathwayScores):
Test cluster-specific genes and pathways via wilcoxon-test across all tissues. For each group j, it is tested feature ~ gr_j vs remaining samples. Combine genes results in loci. Pathways are initially filtered to remove redundant ones and combine Reactome and GO databases. Original genes and pathways are initially corrected for PCs. All tissues can be passed at the same time, script can be parallelized per tissue. Note that the clustering is tissue-specific BUT these scripts test predicted genes and pathways across all tissues.

- *--clusterFile*: complete path to output of clustering script in step 1)
- *--split_tot*: integer indicating the number of groups the *--inputFile* has been split. Used in a large sample size setting e.g. UK Biobank.
- *--type_data_cluster*: type of data used for clustering (tscore, path\_Reactome or path\_GO), corresponds to *--type_data* in step 1)
- *--pvalresFile*: vector, should match *--tissues* entries. TWAS and PALAS results as described in **Input Data**.
- *--ncores*: n. of cores used for parallelization. Parallelizzation per tissue.

#### 2.1.1) Test cluster-specific genes

- *--inputFile*: vector, should match *--tissues* entries. T-score or pathway-scores per tissue as described in **Input Data**.
- *--geneInfoFile*: vector, should match *--tissues* entries. Tab separated file that contains gene location (e.g. resPrior_regEval_allchr.txt from PriLer).
- *--type_data*: type of data that will be tested, default is `tscore` and this is the recommendend mode. However, if path\_Reactome, path\_GO or customPath\_, all pathways in that database will be loaded but without initial filtering and can only be run separately for Reactome, GO or any custom-provided pathway.
- *--pvalcorr_thr*: when combining results in loci for the tested genes, only results with adjusted p-value lower than this threshold will be considered.

```sh
./cluster_associateFeat_corrPCs_run.R \
    --sampleAnnFile \
    --clusterFile \
    --split_tot (default = 0) \
    --inputFile \
    --tissues \
    --type_cluster \
    --functR ./clustering_functions.R \
    --type_data (default "tscore") \
    --type_data_cluster (default "tscore") \
    --type_sim (default "HK") \
    --type_input (default "original") \
    --pvalresFile \
    --geneInfoFile (default NULL) \
    --min_genes_path (default 1) \
    --pval_id (default = 1) \
    --pvalcorr_thr \
    --ncores (default = 5) \
    --outFold 
```
The output includes (saved in *--outFold*):
- **type_data**Original\_corrPCs\_**type_data_cluster**Cluster**type_cluster**\_featAssociation.RData object composed of:
    - inputData: list of loaded *--inputFile*, one per tissue.
    - scaleData: as inputData but scaled per feature and corrected for PCs.
    - res_pval: list of loaded *--pvalresFile*.
    - cl: data frame with clustering partition.
    - tissues: tissues name, entry match inputData and scaleData.
    - covDat: covariates extracted from sampleAnnFile and tested for cluster-specific differences.
    - test_cov: chisq-test or wilcoxon-test for covariates in covDat.
    - test_feat: list of results, one per tissue. Tested each feature in scaleData.
- **type_data**\_corrPCs\_**type_input**\_cluster**type_cluster**\_summary\_geneLoci\_allTissues.txt and **type_data**\_corrPCs\_**type_input**\_cluster**type_cluster**\_summary\_geneLoci\_tissueSpec.txt tab separated tables with results summarized per loci combining all tissues or tissue specific, respectively.

#### 2.1.2) Test cluster-specific pathways
For each tissue, filter pathways based on genes overlap combining both Reactome and GO. It gives priority to pathways with higher coverage and number of genes. It is needed to provide a restricted list of pathways without redundant information.

- *--pvalresFile*: TWAS and PALAS results in .RData object. Used to extract pathway structure.
- *--thr_js*: threshold for jaccard similarity between pathways based on genes, used to clump pathways.

```sh
./filter_pathway_jaccard_sim_run.R \
    --pvalresFile \
    --thr_js (default = 0.2)
    --outFold
```
The output includes (saved in *--outFold*):
-  selected\_pathways\_JSthr**thr_js**.txt: tab separated file. Contains the pathway names to be tested via the next script.


For each tissue, test for cluster-specific pathways. Reactome and GO are merged together. 
- *--inputFold*: vector of complete paths to folder containing pathway-scores. It should match with *--tissues* and should also include selected\_pathways\_JSthr**thr_js**.txt file.
- *--thr_js*: same value used in the previous script.

```sh
./cluster_associatePath_corrPCs_run.R \
    --sampleAnnFile \
    --clusterFile \
    --inputFold \
    --tissues \
    --type_cluster \
    --functR ./clustering_functions.R \
    --type_data_cluster (default "tscore") \
    --type_sim (default "HK") \
    --type_input (default "original") \
    --pvalresFile \
    --pval_id (default 1) \
    --ncores \
    --thr_js (default = 0.2)
    --outFold
```
The output includes (saved in *--outFold*):
- pathOriginal\_filtJS**thr_js**\_corrPCs\_**type_data_cluster**Cluster**type_cluster**\_featAssociation.RData object. Same structure as output of 2.1.1). 

### 2.2) Associate clusters with endophenotypes
Associate the clustering structure with a registered phenotypes on samples. Uses generalized lienar model based on the penotype nature (dependent variables) and corrects for provided covariates. Test 2 models: gr\_i vs remaning samples OR gr\_i vs gr\_j (pairwise) dummy variables as independent. 
- *--sampleAnnFile* MUST include the covariates to correct for.
- *--clusterFile*: complete path to output of clustering script in step 1)
- *--type_data*: NOTE: can also be PCs
- *--type_input*: needed to save output, if refers to version correcting for PCs, add corrPCs\_ in front (e.g. `corrPCs_zscaled`).
- *--risk_score*: if TRUE, the provided phenotypes are predicted gene risk-scores across samples.
- *--rescale_pheno*: if TRUE, continuous phenotype is rescaled.

```sh
./cluster_associatePhenoGLM_run.R \
    --phenoDatFile \
    --phenoDescFile \
    --sampleAnnFile \
    --clusterFile \
    --type_cluster \
    --functR ./clustering_functions.R \
    --type_data \
    --type_sim (default "HK") \
    --type_input (default "original") \
    --risk_score (default FALSE) \
    --rescale_pheno (default FALSE) \
    --outFold
```

The output includes (saved in *--outFold*):
- **type_data**\_**type_input**\_cluster**type_cluster**\_PGmethod\_**type_sim**metric\_phenoAssociation\_GLMpairwise.RData (tests gr\_i vs gr\_j)
- **type_data**\_**type_input**\_cluster**type_cluster**\_PGmethod\_**type_sim**metric\_phenoAssociation\_GLM.RData (tests gr\_i vs remaning samples).
Both .RData objects contain:
    - phenoDat: input --phenoDatFile
    - phenoInfo: input --phenoDescFile
    - cl: data frame with clustering partition
    - covDat: data frame with covariates used
    - bin_reg: summary statistics referring to regression coefficient for group independent variable. 
- **type_data**\_**type_input**\_cluster**type_cluster**_PGmethod_**type_sim**metric_phenoAssociation_GLMpairwise.txt and **type_data**\_**type_input**\_cluster**type_cluster**_PGmethod_**type_sim**metric_phenoAssociation_GLM.txt tab separated file containing `$bin_reg` of the .RData object

### Optional: combine results of multiple phenotype association returns
If mutliple runs of the previous script have been performed to correct for different covariates, the results are combined together in a unique tab separated file. Can also perform a forest plot.
- *--endopFile* vector of .RData oject filed from previous script to be loaded

```sh
./plot_endophenotype_grVSall_run.R \
    --type_cluster_data \
    --type_cluster \
    --type_input \
    --endopFile
    --outFold \
    --forest_plot (default F)\
    --pval_pheno \
    --colorFile 
```
The output includes (saved in *--outFold*):
- **type_cluster_data**\_**type_input**\_cluster**type_cluster**\_PGmethod\_HKmetric\_phenoAssociation\_GLM\_combined.txt

### 2.3) Find differential treatment response among clusters
Test cluster-specific treatment response. It tests the response of phenotypes in `--phenoDatFile` to medication status in `--covDatFile` separately for each group and compares the estimates between each pair of group.
- *--covDatFile*: covariates including treatments in binary format to be tested. The other variables such as PCs are used to corret the model. The model corrects for all the treatments provided.
- *--phenoDescCovFile*: same as --phenoDescFile but referring to --covDatFile. 

```sh
./cluster_treatmentResponseAnalysis_run.R \
    --phenoDatFile \
    --phenoDescFile \
    --phenoDescCovFile \
    --covDatFile \
    --clusterFile \
    --type_cluster \
    --functR ./clustering_functions.R \
    --type_data \
    --type_sim (default "HK") \
    --type_input (default "original") \
    --outFold
```
The output includes (saved in *--outFold*):
- **type_data**\_**type_input**\_cluster**type_cluster**\_TreatResponse\_pairwise.txt summary table for each pair of groups, phenotype and treatment tested. The columns `z_diff` and `pvalue_diff` refers to pairwise difference Z-statistic and p-value.

### 2.4) Drug repositiong based on cluster-specific pathways
Perform drug reporposing based on gene2drug tool implemented in [gep2pep](https://bioconductor.org/packages/release/bioc/html/gep2pep.html) R Bioconductor package. The script searchs for drugs inhibiting group-specific up-regulated pathways or activitating group-specific down-regulated pathways.
- *--pathCluster_file*: output of cluster-specific pathways (2.1.1 script).
- *--atc_file*: csv file obtained from [atcd github](https://github.com/fabkury/atcd) 
- *--cmap_fold*: folder including Cmap data, downloaded from http://dsea.tigem.it/data/Cmap_MSigDB_v6.1_PEPs.tar.gz

```sh
./pathSEA_path_group_run.R \
    --pathCluster_file \
    --atc_file \
    --cmap_fold \
    --type_cluster \
    --outFold
```
The output includes (saved in *--outFold*):
- pathSEA\_corrPCs\_tscoreCluster**type\_cluster**\_featAssociation.txt summary table for all groups and compounds. Includes atc code for the tested compounds, column "type" indicates if up-regulated or donw-regulated pathways in a group were considered.


### 3.1) Project on external cohorts
Project clustering structure into an external cohort not used to derive the original clustering (model cluster). Gene T-scores are pre-processed similar as before and corrected for PCs.
- *--sampleAnnNew_file*: sample file of the new cohort, including covariates
- *--inputFile*: predicted T-scores for the new cohort
- *--sampleAnn_file*: sample file of the cohort already clustered, if NULL "$sampleInfo" from the clustering object is used
- *--clustFile*: .RData object output of step 1)

```sh
./cluster_PGmethod_corrPCs_predict_run.R \
    --inputFile 
    --name_cohort 
    --sampleAnn_file (default NULL)
    --sampleAnnNew_file 
    --clustFile 
    --functR ./clustering_functions.R \
    --type_data (default "tscore") \
    --type_input (default "original") \
    --split_tot (default 0) \
    --type_cluster \ 
    --outFold
```
The output includes (saved in *--outFold*):
- **type\_data**\_corrPCs\_**type\_input**\_predictCluster**type\_cluster**\_PGmethod\_HKmetric.RData R object with the following structure:
    - probability: for each group, indicates the probability of a sample in the external cohort to be projected in that group
    - tot_W_sNN: sparse matrix of shared NN including both the samples in the model clustering and the samples in the external cohort
    - sampleAnn: sampleAnn file of the external cohort
    - data_new: standardized input features used for clustering in the external cohort
    - cl_new: projected clustering structure (assigned as the group having maximum probabiity)
    - res_pval: TWAS assocation for the considered features used for clustering
    - gr_input: mean, sd and coefficient of variation (cv) of each feature across clusters.

### 3.2) Evaluation of projected clustering
Evaluate cluster projection on external cohort compared to model clustering: check percentage of samples in each group, concordance of cluster-specific genes, n.of loci reproduced and test additional endophenotypes (if available). 
- *--model_name*: name of the model clustering cohort(s)
- *--featRel_model*: .RData object output of 2.1.1 for model clustering
- *--clustFile_new*: .RData object output of 3.1 for clustering projection
- *--featRel_predict*: .RData object output of 2.1.1 for clustering projection
- *--geneLoci_summ*: summary\_geneLoci\_allTissues.txt output of 2.1.1 for model clustering

```sh
./cluster_predict_evaluate_run.R \
    --cohort_name \
    --model_name \
    --phenoNew_file (default NULL) \
    --type_cluster \
    --type_data \
    --clustFile \
    --featRel_model  (default NULL) \ 
    --clustFile_new \
    --featRel_predict (default NULL) \
    --functR ./clustering_functions.R \
    --type_input (default "original") \
    --geneLoci_summ (default NULL) \
    --outFold
```
The output includes (saved in *--outFold*):
- **type\_data**\_**type\_input**\_cluster**type\_cluster**\_percentageGropus\_prediction_model**model\_name**.txt percentage of samples in each group in the projected and model clustering 
- **type\_data**\_**type\_input**\_cluster**type\_cluster**\_correlationSpear\_WMWestSign\_Groups\_prediction\_model**model\_name**.txt Spearman correlation between projected and model of wilcoxon-mann-whitney estimates significant in the model clustering.
- **type\_data**\_**type\_input**\_cluster**type\_cluster**\_numberLociRep\_prediction\_model%s.txt Reproducibility of loci significant in the model clustering, uses the most significant gene in loci from the model clustering and compare sign and nominal significance in the projected clustering.
- **type\_data**\_**type\_input**\_cluster**type\_cluster**\_phenoAssociationGLMpairwise\_prediction\_model**model\_name**.RData
and **type\_data**\_**type\_input**\_cluster**type\_cluster**\_phenoAssociationGLM\_prediction\_model**model\_name**.RData same format as in step 2.2

### Optional: Associate projected clusters with endophenotypes
Run endophenotype difference as in 2.2 for a projected clustering. It does not require a phenoInfo file that will be automatically built based on phenotypes nature.
- *--phenoNew_file*: same structure as *--phenoDatFile*
- *--covNew_file*: additional covariates to correct for, if not already available in sample info stored in *--clustFile_new*.

```sh
./cluster_predict_associatePhenoGLM_run.R
    --cohort_name \
    --phenoNew_file \
    --covNew_file (default NULL) \
    --type_cluster \
    --type_data \
    --model_name \
    --clustFile_new \
    --functR ./clustering_functions.R \
    --type_input (default "original") \
    --outFold
```
The output includes (saved in *--outFold*):
- **type\_data**\_**type\_input**\_cluster**type\_cluster**\_phenoAssociationGLMall\_prediction\_model**model\_name**.RData same format as 2.2, the object include both pairwise (gri vs grj) and total (gri vs everything else) comparisons. 

### 4) Associate cluster with gene-risk scores
Predict gene-risk score on the population of interest and detect differences across groups. Gene-risk scores are similar to polygenic risk-score but imputed gene expression and TWAS summary statistics are considered instead of individual genotype and GWAS results. Note that the following script can be lso used for individual pathway-scores and PALAS summary statisitcs but gene usage is set as default.

#### 4.1) Compute features (genes/pathways) correlation
Compute features (genes or pathways) correlation to be used for the initial filtering of features in computing risk-scores. 
A completely different set of samples can be used for this purpose.
- *--inputFile* pathway-score or gene T-scores (.RData) from which the correlation is computed. If pathway-scores, it must be a vector of 2 .RData, one for Reactome and one for GO
- *--sampleAnnFile* samples matching *--inputFile* or a subset used to compute features correlation
- *--split_tot* depends on the dimensionality of input, if 0 a single matrix is loaded

```sh
./correlation_features_run.R \
	--inputFile \
	--sampleAnnFile \
	--tissue_name \
	--split_tot (default 0) \
	--type_data \
	--outFold
```
The output includes (saved in *--outFold*):
- correlation_estimate_<type_data>.RData, R object with correlation matrix and samples info 


#### 4.2) Predict gene risk-scores
- *--genes_to_filter*: .txt files produced by compare_geneExp_matchedDataset_run.R (see below), filter out genes when not harmonized data-sets are considered for TWAS estimates and imputed gene expression
- *--cases_only*: if TRUE, cosider only cases in sampleAnn_file (Dx=1)
- *--scale_rs*: if TRUE risk-scores are standardized
- *--pheno_class_name*: macro name for each RData object including TWAS results in pvalresFile
- *--corrFile*: output of step 4.1)
- *--sqcorr_thr*: squared correlation threshold from --corrFile to clump features based on PriLer R^2.

```sh
./compute_risk_score_corrPCs_run.R
    --genes_to_filter (default NULL) \
    --sampleAnn_file \
    --inputFile \
    --split_tot \
    --n_max (default 50000) \
    --type_cluster \
    --type_data (default "tscore") \
    --cases_only (default FALSE) \
    --functR ./clustering_functions.R \
    --scale_rs (default FALSE) \
    --pheno_class_name \
    --pvalresFile \
    --sqcorr_thr (default 1) \
    --corrFile \
    --min_genes_path (default 1)\
    --outFold
```
The output includes (saved in *--outFold*):
- **type\_data**\_corr2Thr**sqcorr\_thr**\_risk\_score\_relatedPhenotypes.txt, table with first column Individual\_ID and following once predicted gene risk-score for each phenotype having TWAS summary statistic in --pvalresFile.
- **type\_data**\_features\_risk\_score\_corr2Thr**sqcorr\_thr**.txt', table with features (genes) names used to compute gene risk-scores.

#### 4.3) Find cluster-specific differences in gene-risk scores
Application of cluster_associatePhenoGLM_run.R with --risk_score TRUE and --rescale_pheno TRUE

If gene-risk score is predicted among multiple cohorts (with PCs info calculated separately, e.g. CARDIoGRAM), association of gene-RS with each cohort cluster and final meta-analysis computed with the following R script.
- *--name_cohorts*: vector with name of considered cohorts
- *--phenoDatFile*: .txt files including phenotypes to be tested, one per cohort
- *--sampleAnnFile*: .txt files with samples and covariates info, one per cohort
- *--clusterFile*: .RData clustering output, either from projection or actual clustering. Can be a unique file (if all cohorts clustered together) or multiple .RData files one per cohort (cohorts clustered separately)
- *--risk_score*: logical that must be set to TRUE, indicates that the phenotype input is computed as risk score

```sh
./cluster_associatePhenoGLM_multipleCohorts_metaAnalysis_run.R \
    --name_cohorts \
    --phenoDatFile \
    --phenoDescFile \
    --sampleAnnFile \
    --type_cluster \
    --functR ./clustering_functions.R \
    --type_data \
    --type_sim (default "HK") \
    --type_input (default "original") \
    --clusterFile \
    --risk_score (default F) \
    --outFold
```
The output includes (saved in *--outFold*):
- **type_data**\_**type_input**\_cluster**type_cluster**\_PGmethod\_**type_sim**metric\_phenoAssociation\_GLMpairwise_metaAnalysis.RData (tests gr\_i vs gr\_j)
- **type_data**\_**type_input**\_cluster**type_cluster**\_PGmethod\_**type_sim**metric\_phenoAssociation\_GLM_metaAnalysis.RData (tests gr\_i vs remaning samples).
Both .RData objects contain:
    - phenoDat: input --phenoDatFile
    - phenoInfo: input --phenoDescFile
    - cl: data frame with clustering partition
    - covDat: data frame with covariates used
    - single_cohorts: list with statistics referring to regression coefficient for group independent variable, one element per cohort.
    - meta_analysis: meta-analisys results of summary statistics referring to regression coefficient for group independent variable
    
- **type_data**\_**type_input**\_cluster**type_cluster**_PGmethod_**type_sim**metric_phenoAssociation_GLMpairwise_metaAnalysis.txt and **type_data**\_**type_input**\_cluster**type_cluster**_PGmethod_**type_sim**metric_phenoAssociation_GLM_metaAnalysis.txt tab separated file containing `meta_analysis` of the .RData object

#### Optional 4.4): Evaluate gene risk-score prediction
If gene-risk scores are predicted on the same data set that measured the corresponding phenotypes and from which TWAS and PALAS were estimated (e.g. UKBB), it is possible to evaluate the performance of gene-risk scores in approximating the actual phenotypes. The following script computes R2 and F-statistic estiamtes from comparing nested linear models pheno\~risk_score\+covariates and pheno\~risk_score.
- *--riskScore_file*: output of step 4.2)
- *--pheno_file*: phenotype txt file used to compute TWAS and PALAS

```sh
./evaluate_risk_score_run.R \
    --riskScore_file \
    --sampleAnn_file \
    --pheno_file \
    --phenoAnn_file \
    --names_pheno \
    --outFold 
```
The output includes (saved in *--outFold*):
- R2_risk_score_phenotype.txt: summary file with pheno_id, R2 and F-statistic.

F-statistic and R2 across tissues are plot with
```sh
./plot_evaluate_risk_score_run.R \
    --riskScore_eval_file \
    --color_tissues_file \
    --tissues \
    --outFold
```
- *--riskScore_eval_file*: multiple arguments as output of the previous script across tissues, same length as *--tissues*
- *--color_tissues_file*: file with color code per tissue, available in refData/color_tissues.txt 

#### Optional 4.5): Compare results actual endophenotype and gene risk-score cluster-specific differences
Comparison of cluster-specific results from measured endophenotype and gene risk-score. For each phenotype in common, computation of cluter-reliable measure (CRM) given by pheno `F-stat * |beta|` with `beta` being the GLM cluster-specific coefficient from gene risk-score associations.
- *--riskScore_analysis_file*: output of step 4.3) in the form of .RData object. If multiple objects are passed due to multiple run with different endophenotype tables, results are combined together.
- *--endopheno_analysis_file*: output of step 2.2) in .txt format
- *--pval_FDR_pheno*, *--pval_pheno_show*, *--thr_plot*, *--measureGoodness_thr* parameters for plot
- *--color_pheno_file*: color per UKBB phenotype, file available in refData/color_pheno_type_UKBB.txt
- *--R2_pheno_rs_file*: summary table output of step 4.2)
```sh
./compare_endophenotypeAnalysis_clusterRiskScore_run.R \
    --riskScore_analysis_file \
    --endopheno_analysis_file \
    --pval_FDR_pheno (default 0.05) \
    --pval_pheno_show (default 0.001) \
    --thr_plot (default 1e-10) \
    --color_pheno_file \
    --pheno_name \
    --R2_pheno_rs_file \
    --meta_analysis (default F) \
    --measureGoodness_thr (default 3000) \
    --outFold 
```
The output includes (saved in *--outFold*):
- riskScores\_clusterCases\_**pheno_name**\_group\_relatedPheno\_measureGoodnessPred.txt
cluster-specific gene-RS results annotated with F-statistic of gene-RS, actual cluster-specific endophenotypes and cluster-reliable measure.


These files (across tissues) are then used to compute precision across CRM values using beta sign concordance via
```sh
./plot_precision_risk_score_groupSpec_run.R \
    --riskScore_comp_file \
    --color_tissues_file \
    --tissues \
    --outFold
```
The output includes (saved in *--outFold*):
- riskScores\_clusterCases\_group\_relatedPheno\_npheno\_and\_precision.txt
across CRM threshold gives info on n. of phenotypes passing CRM and having the same sign for beta GLM association and resulting precision.
***

## Multiple cohorts (tissue-specific)

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
- **tissue_name**\_filter\_genes\_matched\_datasets.txt, the column "keep" is a logical indicating whether a gene passed the correlation threshold and is reliable in both models.

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
- **tissue_name**\_filter\_**path_type**\_path\_matched\_datasets.txt, includes only pathways in common and the column "keep" is a logical indicating whether a pathway passed the correlation threshold.

### 0) Detect outliers
Prior to clustering computation combining all cohorts, detect samples outliers based on UMAP representation. Samples more distant than 6 standard deviation from the median are excluded (union considering the first two dimension). The script is run separately for each tissue.
- *--name_cohorts*: vector with name of considered cohorts
- *--sampleAnnFile*: .txt files with samples and covariates info, one per cohort
- *--genes_to_filter*: .txt output of previous optional step. If not NULL, only genes in "keep" are considered for the analysis. 
- *--pval_id*: integer indicating the index of the phenotype to be considered in *--pvalresFile* that can include more than one phenotype.
- *--corr_thr*: correlation threshold to clump features. If 1 include all the features.
- *--min_genes_path*: minimum number of genes forming a pathway. Used to filter pathways for clustering with *--type_data* is path_Reactome or path_GO.

```sh
./detect_outliers_corrPCs_multipleCohorts_run.R \
    --inputFile \
    --name_cohorts \
    --sampleAnnFile \
    --geneRegionFile (default NULL) \
    --tissues_name \
    --genes_to_filter (default NULL) \
    --exclude_MHC (default FALSE) \
    --type_cluster \
    --corr_thr (default -1) \
    --split_tot (default 0) \
    --pvalresFile \
    --pval_id (default 1) \
    --functR ./clustering_functions.R \
    --type_data \
    --type_sim (default "HK") \
    --type_input (default "original") \
    --min_genes_path (default 1) \
    --outFold
```
The output includes (saved in *--outFold*):
- **type_data**\_corrPCs\_**type_input**\_cluster**type_cluster**\_PGmethod\_umap\_oultiers.txt containing the sample ids of outliers to exclude. If no outliers are detected, no file is saved.

Combine list of outliers across multiple tissues and/or clustering versions (e.g. different corr_thr).
- *--sampleFiles*: vector of .txt files with outlier samples

```sh
./combine_outliers_cluster_run.R \
    --sampleFiles \
    --type_cluster \
    --type_data \
    --type_input (default "original") \
    --outFold
```
The output includes (saved in *--outFold*):
- samples\_to\_remove\_outliersUMAP\_**type_data**\_**type_input**\_cluster**type_cluster**.txt includes all samples considered as outliers across multiple tissues/clustering versions.

### 1) Clustering
Cluster individuals based on tissue-specific imputed expression. The individuals are split across K cohorts. These cohorts are concatenated, genes are clumped based on the correlation computed from the concatenated dataset. Each gene is standardized, corrected for PCs and rescaled by TWAS Z-statistic. Note that correcting for PCs after the datasets have been concatenated is reasonable only if PCs are computed merging the genetic imput of the multiple cohorts. 
- *--inputFile*: vector of txt files (predictedTscores.txt) or common name of split .RData files (predictedTscore_splitGenes), one per cohort
- *--sampleAnnFile*: .txt files with samples and covariates info, one per cohort
- *--sampleOutFile*: .txt tab separated file, list of samples to be removed from clustering, columns must include Individual_ID. It can be the output of step 0)

```sh
./cluster_PGmethod_corrPCs_multipleCohorts_run.R \
	--inputFile \
	--genes_to_filter (default NULL) \
	--name_cohorts \
	--sampleAnnFile \
	--sampleOutFile \
	--geneRegionFile (default NULL) \
	--tissues_name \
	--exclude_MHC (default FALSE) \
	--type_cluster \
	--split_tot (default 0) \
	--pvalresFile \
	--pval_id (default 1) \
	--corr_thr (default -1) \
	--functR ./clustering_functions.R \
	--type_data (default tscore) \
	--type_sim (default HK) \
	--type_input (default original)
	--kNN_par (default 20) \
    --cluster_method (default "leiden") \
	--min_genes_path (default 1) \
	--outFold
```
The output includes (saved in *--outFold*):
- **type_data**\_corrPCs\_**type_input**\_cluster**type_cluster**\_PGmethod\_**type_sim**metric.RData object with the same structure as step 1) in "Single cohort". In addition it includes 
	- sampleInfo: complete sample annotation with cohort division
	- sampleOutliers: IDs of removed samples
	- umap: First 2 UMAP coordinates
- **type_data**\_corrPCs\_**type_input**\_cluster**type_cluster**\_PGmethod\_**type_sim**metric_minimal.RData similar object but only including
	- cl_best
	- samples_id
	- feat 
	- test_cov

### 2.1) Associate clusters with molecular features (genes/pathwayScores):
Test cluster-specific genes and pathways via wilcoxon-test across all tissues. Similar to Step 2.1) of "Single cohort" but associations tested concatenating all cohorts. Differences are in 
- *--sampleAnnFile*: vector of .txt files refering to sample annotation, one per cohort
- *--inputFold*: vector of input (tscore or pathway-score) fold, one per tissue (common part across cohorts)
- *--name_cohorts*: vector of cohort names

#### 2.1.1) Test cluster-specific genes
Similar to 2.1.1) of "Single cohort". Cna also test differences at the level pf pathways, provided that the .txt input file is given. Pathways are NOT filtered and all the given ones are tested. Differences are in
- *--additional_name_file*: common final part referring to tscore or pathway-score (e.g. predictedTscores.txt) common across cohorts. Note that the input file complete name is built as `name_file <- paste0(inputFold[id_t], name_cohorts[c_id], additional_name_file)`

```sh
./cluster_associateFeat_corrPCs_multipleCohorts_run.R \
	--sampleAnnFile \
	--clusterFile \
	--split_tot (default 0) \
	--name_cohorts \
	--inputFold \
	--additional_name_file \
	--tissues \
	--type_cluster \
	--functR ./clustering_functions.R \
	--type_data (default "tscore") \
	--type_data_cluster (default "tscore") \
	--type_sim (default "HK") \
	--type_input (default "original") \
	--pvalresFile \
	--geneInfoFile \
	--min_genes_path (default 1) \
	--pval_id (default 1) \
	--pvalcorr_thr (default 0.05)
	--ncores (default 5) \
	--outFold
```
The output includes (saved in *--outFold*):
- **type_data**Original\_corrPCs\_**type_data_cluster**Cluster**type_cluster**\_featAssociation.RData object composed of:
    - inputData: list of loaded *--inputFile*, one per tissue.
    - scaleData: as inputData but scaled per feature and corrected for PCs.
    - res_pval: list of loaded *--pvalresFile*.
    - cl: data frame with clustering partition.
    - tissues: tissues name, entry match inputData and scaleData.
    - covDat: covariates extracted from sampleAnnFile and tested for cluster-specific differences.
    - test_cov: chisq-test or wilcoxon-test for covariates in covDat.
    - test_feat: list of results, one per tissue. Tested each feature in scaleData.
- **type_data**\_corrPCs\_**type_input**\_cluster**type_cluster**\_summary\_geneLoci\_allTissues.txt and **type_data**\_corrPCs\_**type_input**\_cluster**type_cluster**\_summary\_geneLoci\_tissueSpec.txt tab separated tables with results summarized per loci combining all tissues or tissue specific, respectively.

#### 2.1.2) Test cluster-specific pathways
For each tissue, filter pathways based on genes overlap combining both Reactome and GO. It gives priority to pathways with higher coverage and number of genes. It is needed to provide a restricted list of pathways without redundant information. Same as step 2.1.2 in "Single cohort"
```sh
./filter_pathway_jaccard_sim_run.R \
    --pvalresFile \
    --thr_js (default = 0.2)
    --outFold
```
The output includes (saved in *--outFold*):
-  selected\_pathways\_JSthr**thr_js**.txt: tab separated file. Contains the pathway names to be tested via the next script.

For each tissue, test for cluster-specific pathways. Reactome and GO are merged together. Same as step 2.1.2) of "Single cohorts" but test concatanating all cohorts. Differences are in
- *--additional_name_file*: common final part referring pathway-score location including both GO and Reactome pathways (e.g. /devgeno0.01_testdevgeno0/) and common across cohorts. Note that the input file complete name is built as `name_file <- paste0(inputFold[id_t], name_cohorts[c_id], additional_name_file)` and `inputFile = sprintf('%s/Pathway_Reactome_scores.txt', name_file)`

```sh
./cluster_associatePath_corrPCs_multipleCohort_run.R \
	--sampleAnnFile \
	--clusterFile \
	--name_cohorts \
	--inputFold \
	--additional_name_file \
	--tissues \
	--type_cluster \
	--functR ./clustering_functions.R \
	--type_data \
	--type_data_cluster \
	--type_sim (default "HK") \
	--type_input (default "original") \
	--pvalresFile \
	--pval_id (default 1) \
	--ncores (default 5) \
	--thr_js (default 0.2) \
	--path_filt_file \
	--outFold
```
The output includes (saved in --outFold):
- pathOriginal\_filtJS**thr_js**\_corrPCs\_**type_data_cluster**Cluster**type_cluster**\_featAssociation.RData object. Same structure as output of 2.1.1). 

#### 2.4) Drug repositiong based on cluster-specific pathways
Perform drug reporposing, same as in step 2.4) of "Single cohort"

### 3) Project on external cohorts
Project clustering structure into an external cohort not used to derive the original clustering (model cluster). Scripts in steps 3) for "Single cohort" can be applied in the same way.

### 4) Associate cluster with gene-risk scores
Predict gene-risk score on the population of interest and detect differences across groups (see "Single cohort" step 4) for more details).

#### 4.1) Compute features (genes/pathways) correlation
Same script as in "Single cohort" step 4.1), use external single cohort samples for this purpuse.

#### 4.2) Predict gene risk-scores
Similar to step 4.2) in "Single cohort". It computes gene risk scores across all samples in multiple cohorts, after having clumped genes and corrected the remaining genes for PCs.
- *--inputFile*: vector including file direction or complete path for gene T-scores, same length as --name_cohorts
```sh
./compute_risk_score_corrPCs_multipleCohorts_run.R \
	--genes_to_filter (default NULL) \
	--sampleAnn_file \
	--name_cohorts \
	--inputFile \
	--split_tot (deafult 0)
	--functR ./clustering_functions.R \
	--type_data \
	--cases_only (default FALSE) \
	--scale_rs (default FALSE) \
	--pheno_class_name \
	--pvalresFile \
	--sqcorr_thr (default 1) \
	--n_max (default 50000) \
	--corrFile \
	--min_genes_path (default 1) \
	--outFold
```
The output includes (saved in *--outFold*):
- **type\_data**\_corr2Thr**sqcorr\_thr**\_risk\_score\_relatedPhenotypes.txt, table with first column Individual\_ID and following once predicted gene risk-score for each phenotype having TWAS summary statistic in --pvalresFile.
- **type\_data**\_features\_risk\_score\_corr2Thr**sqcorr\_thr**.txt', table with features (genes) names used to compute gene risk-scores.

#### 4.3) Find cluster-specific differences in gene-risk scores
Application of cluster_associatePhenoGLM_run.R with --risk_score TRUE and --rescale_pheno TRUE





