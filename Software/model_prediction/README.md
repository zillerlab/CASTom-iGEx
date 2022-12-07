# Imputation of gene expression, computation of individual pathway-scores and association with trait 
CASTom-iGEx (Module 2) is a pipeline (R based) that uses trained model to 
1. predicts gene expression from genotype-only datasets and convert them to T-scores and pathway scores (two versions are available dependening of data dimensionality); 
2. perform association analsysis of genes (TWAS) and pathways (PALAS) with the traits of interest.

## Input Files
- **Genotype matrix** (*--genoDat_file*): dosages for each chromosome (compressed txt) without variants name/position (variants x samples).  *NOTE: SNPs must match (or be a subset of) the train genotype data, file must end with chr<>_matrix.txt.gz*
- **Phenotype matrix** (*--phenoDat_file*): columns must contain `Individual_ID` plus any phenotype to test the association (phenotypes + 1 x samples). This matrix can include multiple phenotypes to be tested. 
- **Phenotype description (*--phenoAnn_file*)**: csv file, rows refers to phenotypes to be tested. Columns must include: `pheno_id`, `FieldID`, `Field`,  `transformed_type`;  `pheno_id` is used to match columns name in Phenotype matrix, `transformed_type` is a charachter defining the type of data. Inspired by PHESANT for UKBiobank, possible values of `transformed_type` are 
    - "CONTINUOUS" (gaussian regression)
    - "CAT_SINGLE_UNORDERED", "CAT_SINGLE_BINARY", "CAT_MUL_BINARY_VAR" for binary (binomial regression)
    - "CAT_ORD" for ordinal (ordered logistic regression)
- **Covariate matrix** (*--covDat_file*): covariates to correct for in the association analysis (covariats + IDs x samples). Columns must contain `Individual_ID` and `genoSample_ID` to match genotype plus covariates to correct for in the phenotype association. Column `Dx` (0 control 1 case) is optional, if present is used to build the reference set when computing T-scores. *Note: samples in genotype and phenotype matrix are matched based on covariate matrix*
- **Reactome Pathway annotation** (*--reactome_file*): .gmt file can be downloaded from https://reactome.org/download-data/ (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **GO Pathway annotation** (*--GOterms_file*): .RData file, can be obtained using *Annotate_GOterm_run.R*, each pathway is a entry in the list with `GOID` `Term` `Ontology` `geneIds` elemnets (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **Custom pathway** (*--pathwayStructure_file*): .RData file, similar to GO structure, each pathway is a list entry with `name` and `geneIds` elements. Available for WikiPathways (2019) in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/)

## Workflow
### 1) Predict gene expression
From previously trained PriLer tissue-specific model (Module 1), predict gene expression based on genotype-only data set

*--InfoFold* is the folder with gene-snp distance matrix ENSEMBL_gene_SNP_2e+5_chr<>_matrix.mtx, *--outTrain_fold* is the folder with overall results from PriLer

```sh
./Priler_predictGeneExp_run.R \
    --genoDat_file \
    --covDat_file \
    --cis_thres (default = 200000) \
    --outFold \
    --outTrain_fold  \
    --InfoFold \
    --no_zip (default F)
```
*NOTE: can be split for subset of samples (depends on covDat_file), genes are NOT filtered*

The output includes (saved in *--outFold*):
- predictedExpression.txt.gz 

### 1') Alternative: Predict gene expression on not harmonized data
Similar as before BUT the genotype-only data is not initially harmonized with the reference panel. This scenario happens if it's necessary to use a previously trained model on a new data, however the best solution would be to initially harmonize the reference panel and the genotype-only data to properly account for the LD structure of the available variants. Prior to this step, genotype-only must be filtered to include only variants in the reference panel with the same REF and ALT annotation. The gene expression is imputed based on variants intersection.

```sh
./Priler_predictGeneExp_smallerVariantSet_run.R \
    --genoDat_file \
    --covDat_file \
    --cis_thres (default = 200000) \
    --outFold \
    --outTrain_fold  \
    --InfoFold \
    --no_zip (default F) \
    --genoInfo_model_file \
    --genoInfo_file \
```

- *--genoInfo_file* and *--genoInfo_model_file* are the files having variant annotation of the new genotype-only dataset and from PriLer model that must include columns POS, REF and ALT. The format file include the complete path excluding the last part that should end with chr<>.txt
- rows in *--genoDat_file* must match rows in *--genoInfo_file*

The output includes (saved in *--outFold*):
- predictedExpression.txt.gz 

***
Based on data dimension, the next scripts are divided in two parts. If sample size <= 10,000 follow "Small dataset" part, otherwise "Large dataset".
***

### 2) Small dataset: T-scores and Pathway-scores computation
Predicted gene expression is converted into T-scores and combined into Pathway-scores. T-scores are computed using a subset of control samples (`Dx == 0`) as reference set if column `Dx` is present in *--covDat_file*, otherwise a random subset of samples is selected. 

*--input_file* is predictedExpression.txt.gz of the previous step, 

```sh
./Tscore_PathScore_diff_run.R \
    --input_file \
    --reactome_file \
    --GOterms_file \
    --originalRNA (default = F) \
    --thr_reliableGenes (default = c(0.01, 0)) \
    --covDat_file \
    --nFolds (default = 20) \
    --outFold
```
The output includes (saved in *--outFold*):
- predictedTscores.txt: gene T-scores (used HUGO nomenclature)
- Pathway_Reactome_scores.txt: pathway scores based on Reactome
- Pathway_GO_scores.txt: pathway scores based on GO

Pathway scores can include pathway made of the same genes but with different names (filtered in the association steps)

### 2') Small dataset: Pathway-scores computation for custom gene list
Based on already computed T-scores, create pathway-scores for a custom .RData object containing gene sets (example WikiPathways)

*--tscore_file* is - predictedTscores.txt: gene T-scores of the previous step, *--geneSetName* custom name for gene sets databse

```sh
./pathScore_customGeneList_run.R \
    --pathwayStruct_file \
    --tscore_file \
    --sampleAnn_file \
    --geneSetName \
    --abs_tscore (default = F) \
    --outFold
```
The output includes (saved in *--outFold*):
- Pathway_<*geneSetName*>_scores.txt

### 3) Small dataset: Association with phenotype of T-score and pathways
T-scores and pathway-scores are tested for association with phenotypes. The regression type depends on the nature of the phenotype (gaussian, binary and ordinal logistic). Redundant pathways composed of the same gene sets are removed keeping the one with lower number of annotated total gene. Genes/pathways are corrected for multiple testing.
- *--thr_reliableGenes* MUST be the same as the filtering criteria previously applied, 
- *--inputFold* contains pathways and t-scores results, 
- *--covDat_file* includes covariates to filter for, 
- *--sampleAnn_file* same sample list of *covDat_file* but can exclude actual covariates, 
- *--names_file* name for the phenotype group to be tested, 
- *--geneAnn_file* resPrior_regEval_allchr.txt from PriLer
- *--cov_corr* if TRUE (default) corrects for the covariates available in covDat_file

```sh
./pheno_association_smallData_run.R \
    --reactome_file \
    --GOterms_file \
    --sampleAnn_file \
    --thr_reliableGenes (default = c(0.01, 0)) \
    --covDat_file \
    --phenoDat_file \
    --names_file \
    --phenoAnn_file \
    --cov_corr (default = T) \
    --ncores (default = 0) \
    --geneAnn_file \
    --functR ./pheno_association_functions.R \
    --outFold
```
The output includes (saved in *--outFold*):
- pval_<*names_file*>_covCorr.RData list composed of 
    - pheno: phenoAnn info
    - tscore (list, each entry refers to a phenotype): summary statistics of association from tscores for each reliable genes
    - pathScore_reactome (list, each entry refers to a phenotype): summary statistics of association from pathscore (Reactome)
    - pathScore_GO (list, each entry refers to a phenotype): summary statistics of association from pathscore (GO)
    - info_pathScore_reactome (list, each entry refers to a phenotype): for each pathway in Reactome, its summary statistics and those of the genes belonging to the pathway
    - info_pathScore_GO (list, each entry refers to a phenotype): for each pathway in GO, its summary statistics and those of the genes belonging to the pathway


### 4) Small dataset: Association with phenotype of custom pathways
Same as before but for custom gene sets. It requires the association between phenotypes and T-scores to be complete (previous step).

```sh
./pheno_association_smallData_customPath_run.R \
    --pathwayStructure_file  \
    --sampleAnn_file \
    --thr_reliableGenes (default = c(0.01, 0)) \
    --covDat_file \
    --phenoDat_file \
    --names_file \
    --phenoAnn_file \
    --cov_corr (default = T) \
    --ncores (default = 0) \
    --geneAnn_file \
    --functR ./pheno_association_functions.R \
    --geneSetName \
    --abs_tscore (default F) \
    --not_rm_samepath (default F) \
    --outFold
```

The output includes (saved in *--outFold*):
- pval_<*names_file*>_covCorr_customPath_<*geneSetName*>.RData (same structure as previous step)

### 5) Small dataset: Meta-analysis for multiple cohorts T-scores and pathways

Combine results from multiple cohorts (harmonized) via meta-analsys (inverse variance weighted method). 
- *--res_cohorts* .RData files from pheno_association_smallData_run.R, one for each cohort
- *--thr_het* threshold for p-value from Cochrane’s Q statistic, if p-value <= thr_het estimates are computed via random-effects

```sh
./pheno_association_metaAnalysis_run.R \
    --res_cohorts
    --name_cohort 
    --phenoDatFile_cohorts
    --cov_corr (defualt T)
    --phenoName
    --thr_het (defualt 0.001)
    --reactome_file
    --GOterms_file
    --outFold
```
The output includes (saved in *--outFold*):
- pval_<*phenoName*>_covCorr.RData (same structure as previous step)
- phenoInfo_<*phenoName*>_cohorts.txt (tab separated file with n. cases/controls for each sample)

### 6) Small dataset: Meta-analysis for multiple cohorts custom pathways 
Same as before but for custom gene sets.

```sh
./pheno_association_customPath_metaAnalysis_run.R \
    --res_cohorts 
    --name_cohort
    --phenoDatFile_cohorts
    --cov_corr (default T)
    --phenoName
    --thr_het (defualt 0.001)
    --pathwayStructure_file
    --geneSetName 
    --abs_tscore (default T)
    --outFold 
```
The output includes (saved in *--outFold*):
- pval_<*phenoName*>_covCorr.RData 
- phenoInfo_<*phenoName*>_cohorts.txt (tab separated file with n. cases/controls for each sample)

***

### 1) Large dataset: preliminary

Predicted gene expression had been executed for split set of samples. For each of them keep only gene such that dev_geno > 0.01 and test_dev_geno > 0.
- *inputFold* folder including resPrior_regEval_allchr.txt from PriLer
- *outFold* folder including predicted gene expression (split per sample groups)
- *split_tot* number of subgroup the samples have been split on

```sh
bash Combine_filteredGeneExpr.sh inputFold outFold split_tot
```

The output includes:
- predictedExpression_filt.txt.gz filtered predicted gene expression for all samples

### 2) Large dataset: T-scores computation
Perform differential gene expression using t-statistic for each sample with respect to subset of samples considered as reference. The reference is usually a subset of control samples. However if column Dx is not present in the covariate mat (covDat_file), a subset of samples is randomly chosen.

**NOTE**: computationally heavy. All samples considered together, process is split across genes (--split_gene_id)

- *input_file* vector containing full path to split gene expression files
- *split_tot* number of subgroup the genes will be split on, each split will produce n_samples x n_genes_split matrices

```sh
./Tscore_splitGenes_run.R \
    --input_file \
    --nFolds (default 40) \
    --perc_comp (default 0.5) \
    --ncores (default 10) \
    --covDat_file \
    --outFold \
    --split_gene_id \
    --split_tot (default 100)
```

The output includes (saved in *--outFold*):
- predictedTscore_splitGenes{i}.RData for each subset of genes, all samples included

### 3) Large dataset: PathScore computation
Combine T-scores into Pathway scores using as annotation Reactome and GO.
**NOTE**: computationally heavy, for ~ 300,000 samples it required --mem-per-cpu=35G and --cpus-per-task=10

- *input_file* common path to .RData object predicted Tscores ({i}.RData part excluded)
- *split_tot* MUST be tha same value used in the previous script

```sh
./PathwayScores_splitGenes_run.R \
    --ncores (default 10) \
    --input_file  \
    --covDat_file \
    --outFold \
    --split_tot (default 100) \
    --reactome_file \
    --GOterms_file \
    --skip_reactome (default F)
```

The output includes (saved in *--outFold*):
- Pathway_Reactome.RData: pathway scores based on Reactome
- Pathway_GO_scores.RData: pathway scores based on GO

### 3') Large dataset: Pathway-scores computation for custom gene list
Combine T-scores into Pathway scores for a custom .RData object containing gene sets (example WikiPathways)
**NOTE**: computationally heavy, for ~ 300,000 samples it required --mem-per-cpu=35G and --cpus-per-task=10

```sh
./PathwayScores_splitGenes_customGeneList_run.R \
    --ncores (default 10) \
    --input_file  \
    --covDat_file \
    --outFold \
    --split_tot (default 100) \
    --pathwayStruct_file \
    --geneSetName
```

The output includes (saved in *--outFold*):
- Pathway_<*geneSetName*>.RData 

### 4) Large dataset: Association
Perform association with a phenotype, it is divided in multiple steps to reduce computational burden

#### 4.1) Large dataset: Association: input preparation
Prepare for association with phenotypes, it creates summary info files for genes and pathways and split pathway .RData in smaller .RData objectes to be tested in parallel. Pathways computed from the same genes are filtered to keep the ones composed of a lower number of genes (original annotation).
- *--inputFold* folder including pathway and T-scores from step 3)
- *--geneAnn_file* resPrior_regEval_allchr.txt from PriLer
- *--skip_tscore_info* boolen, set to TRUE if preparation for T-score has been done, e.g. when applying the script to custom pathway
- *--split_tot* number of subgroup the genes were split on, the same value will be used to split pathways into that amount of groups.

```sh
./pheno_association_prepare_largeData_run.R \
    --inputFold (default NULL) \
    --reactome_file (default NULL) \
    --GOterms_file (default NULL) \
    --pathwayCustom_file (default NULL) \
    --pathwayCustom_name (default NULL) \
    --sampleAnn_file \
    --geneAnn_file \
    --thr_reliableGenes (default = c(0.01, 0)) \
    --split_tot (defualt 100) \
    --skip_tscore_info (default FALSE) \
    --outFold 
```
The output includes (saved in *--outFold*):
- tscore_info.RData, pathScore_Reactome_info.RData, pathScore_GO_info.RData, pathScore_<pathwayCustom_name>_info.RData tables with genes and pathways to be tested
- Pathway_Reactome_scores_splitPath<i>.RData, Pathway_GO_scores_splitPath<i>.RData, Pathway_<pathwayCustom_name>_scores_splitPath<i>.RData matrices of subgroup of pathways to be tested, i goes from 1 to *split_tot*

#### 4.2) Large dataset: Association: test gene T-scores
Test gene T-scores association with phenotype via generalized linear model. The script is run for each group of genes (split **i**) and for each provided phenotype. 
- *--inputFile* complete path to predictedTscore_splitGenes{i}.RData,
- *--inputInfoFile* complete path to tscore_info.RData obtained from the previous script,
- *--split_tot* MUST be the same value used in Tscore_splitGenes_run.R,
- *--split_gene_id* value between 1 and split_tot indicating the group of genes to be tested,
- *--sampleAnn_file* same sample list of *covDat_file* but can exclude actual covariates,
- *--cov_corr* if TRUE (default) corrects for the covariates available in covDat_file,
- *--names_file* name for the phenotype group to be tested.
Note that multiple files can be passed to *phenoDat_file* and *covDat_file* that MUST be of the same lenght of the provided *names_file* vector. Each entry will refer to a group of phenotypes to be tested and the covariates to correct for. 

```sh
./pheno_association_tscore_largeData_run.R \
    --inputFile \
    --inputInfoFile \
    --split_tot (defualt 100) \
    --split_gene_id \
    --covDat_file \
    --sampleAnn_file \
    --phenoDat_file \
    --phenoAnn_file \
    --cov_corr (default TRUE)\
    --names_file \
    --functR ./pheno_association_functions.R \
    --ncores \
    --outFile \
```
The output is:

#### 4.3) Large dataset: Association: test pathway-scores

```sh
./pheno_association_pathscore_largeData_run.R \
```

#### 4.4) Large dataset: Association: combine results

```sh
./pheno_association_combine_largeData_run.R \
```

#### 4.4') Large dataset: Association: combine results from custom pathway database
```sh
./pheno_association_combine_largeData_customPath_run.R \
```

***



