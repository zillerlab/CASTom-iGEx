# Imputation of gene expression, computation of individual pathway-scores and association with trait 
CASTom-iGEx (Module 2) is a command-line tool that uses trained model to predicts gene expression from genotype-only datasets and convert them to T-scores and pathway scores. Two versions are available dependening of data dimensionality. Prediction step contains also scripts to perform association with trait of interest as well as perform mendelian randomization across traits based on based on pathway and genes association.

## Input Files
- **Genotype matrix** (*--genoDat_file*): dosages for each chromosome (compressed txt) without variants name/position (variants x samples).  *NOTE: SNPs must match with the train genotype data, file must end with chr<>_matrix.txt.gz*
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
### Predict gene expression
From previously trained PriLer tissue-specific model (Module 1), predict gene expression based on genotype-only dataset

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
*NOTE: can be splitted for subset of samples (depends on covDat_file), genes are NOT filtered*

The output includes (saved in *--outFold*):
- predictedExpression.txt.gz 

### Alternative: Predict gene expression on not harmonized data
Similar as before BUT the genotype-only data is not initially harmonized with the reference panel. This scenario happens if it's necessary to use a previously trained model on a new data, however the best solution would be to initially harmonize reference pane and genotype-only data so to preoperly account for the LD structure of the available variants. Prior to this step is necessary that the genotype-only data has been filtered to include only variants from the reference panel genotype data (also having same REF and ALT annotation). The gene expression are imputed based on variants intersection.

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

### Small dataset: T-scores and Pathway-scores computation
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

### Small dataset: Pathway-scores computation for custom gene list
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

### Small dataset: Association with phenotype of T-score and pathways
T-scores and pathway-scores are tested for association with phenotypes. The regression type depends on the nature of the phenotype (gaussin, binary and ordinal logistic). Redundant pathways composed of the same gene sets are removed keeping the one with lower number of annotated total gene. Genes/pathways are corrected for multiple testing.

*--thr_reliableGenes* MUST be the same as the filtering criteria previously applied, *--inputFold* contains pathways and t-scores results, *--covDat_file* includes covariates to filter for, *--sampleAnn_file* same sample list of *covDat_file* but can exclude actual covariates, *--names_file* name for the phenotype group to be tested, *--geneAnn_file* resPrior_regEval_allchr.txt from PriLer

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
- pval_<*names_file*>_covCorr.RData/pval_<*names_file*>.RData list composed of 
    - pheno: phenoAnn info
    - tscore (list, each entry refers to a phenotype): summary statistics of association from tscores for each reliable genes
    - pathScore_reactome (list, each entry refers to a phenotype): summary statistics of association from pathscore (Reactome)
    - pathScore_GO (list, each entry refers to a phenotype): summary statistics of association from pathscore (GO)
    - info_pathScore_reactome (list, each entry refers to a phenotype): for each pathway in Reactome, its summary statistics and those of the genes belonging to the pathway
    - info_pathScore_GO (list, each entry refers to a phenotype): for each pathway in GO, its summary statistics and those of the genes belonging to the pathway


### Small dataset: Association with phenotype of custom pathways
Same as before btu for custom gene sets. It requires the association between phenotypes and T-scores to be complete (previous step).

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
- pval_<*names_file*>_covCorr_customPath_<*geneSetName*>.RData / pval_<*names_file*>_customPath_<*geneSetName*>.RData (same as previous step)

### Small dataset: Meta-analysis for multiple cohorts T-scores and pathways
pheno_association_metaAnalysis_run.R

### Small dataset: Meta-analysis for multiple cohorts custom pathways 
pheno_association_customPath_metaAnalysis_run.R

***
***

### 1) Large dataset: preliminary
Predicted gene expression had been executed for split set of samples. For each of them keep only gene such that dev_geno>0.01 and test_dev_geno>0.

#### Usage
> bash Combine_filteredGeneExpr.sh inputFold outFold split_tot

- inputFold: folder including training evaluation model
- outFold: folder including predicted gene expression
- split_tot: number of subgroup the samples have been split on

The output includes:
- split{i}_predictedExpression_filt.txt filtered predicted gene expression for each sample subset
- predictedExpression_filt.txt.gz filtered predicted gene expression for all samples

### 2) Large dataset: T-scores computation
Perform differential gene expression using t-statistic for each sample with respect to subset of samples considered as reference. The reference is usually a subset of control samples, however if column Dx is not present in the covariate Matrix file a subset of samples is randomly chosen.
*NOTE: computationally heavy. All samples considered together, process is split across genes*
#### Usage
> Rscript Tscore_splitGenes_run.R --input_file --nFolds (default 20) --perc_comp (default 0.5) --ncores (default 10) --covDat_file --outFold --split_gene_id --split_tot (default 100)

- *input_file*: vector containing full path to split gene expression files

The output includes: 
- predictedTscore_splitGenes{i}.RData for each subset of genes, all samples included

### 3) Large dataset: PathScore computation
Combine T-scores into Pathway scores using as annotation Reactome and GO.
*NOTE: computationally heavy*
#### Usage
> Rscript PathwayScores_splitGenes_run.R/PathwayScores_splitGenes_customGeneList_run.R --ncores (default 10) --input_file  --covDat_file  --outFold --split_tot (default 100) --reactome_file --GOterms_file --skip_reactome (default F)
- *input_file*: common path to .RData object predicted Tscores ({i}.RData part excluded)

The output includes:
- Pathway_Reactome/GO_scores.RData: pathway scores for each pathways and samples
- Pathway_Reactome/GO_pvalues.RData: pathway pvalues (t.test) for each pathways and samples
- PathwaySummaryPvalues_cases_Reactome/GO.txt: Fisher’s combined probability pvalue consider only cases
- PathwaySummaryPvalues_cases_Reactome/GO.txt: Fisher’s combined probability pvalue consider only controls

#### 4) Large dataset: Association
- pheno_association_prepare_largeData_run.R
- pheno_association_tscore_largeData_run.R
- pheno_association_pathscore_largeData_run.R
- pheno_association_combine_largeData_run.R

***
***

### Initial filtering if datasets are not harmonized 

compare_geneExp_matchedDataset_run.R
compare_pathScore_matchedDataset_run.R

### Correlation based on gene and pathway association of a trait of interest with multiple endophenotypes 

correlation_pheno_relatedPheno_run.R

### Mendelian randomization based on genes and pathways association 

correlation_features_run.R

#### (direct)
mendelianRand_pheno_relatedPheno_run.R

#### (reverse)

mendelianRand_reverse_pheno_relatedPheno_run.R


