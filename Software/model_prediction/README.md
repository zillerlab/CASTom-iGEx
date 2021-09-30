# Imputation of gene expression, computation of individual pathway-scores and association with trait 
CASTom-iGEx (Module 2) is a command-line tool that uses trained model to predicts gene expression from genotype-only datasets and convert them to T-scores and pathway scores. Two versions are available dependening of data dimensionality. Prediction step contains also scripts to perform association with trait of interest as well as perform mendelian randomization across traits based on based on pathway and genes association.

## Requirements
To run the prediction the following R packages are required:
 - argparse
 - Matrix
 - bigmemory
 - pryr
 - qvalue
 - matrixStats
 - parallel
 - doParallel
 - data.table
 - PGSEA
 - limma
 - nnet
 - lattice
 - MASS
 - lmtest


## Input Files
- **Genotype matrix** (*--genoDat_file*): dosages for each chromosome (compressed txt) without variants name/position (variants x samples).  *NOTE: SNPs must match with the train genotype data, file must end with chr<>_matrix.txt.gz*
- **Phenotype matrix**: columns must contain `Individual_ID` plus any phenotype to test the association (phenotypes + 1 x samples). This matrix can include multiple phenotypes to be tested. 
- **Phenotype description**: rows refers to phenotypes to be tested. Columns must include: `pheno_id`, `FieldID`, `Field`,  `transformed_type`;  `pheno_id` is used to match columns name in Phenotype matrix, transformed_type is a charachter defining the type of data (continous, binary ecc.)
- **Covariate matrix** (*--covDat_file*): covariates to correct for in the association analysis (covariats + IDs x samples). Columns must contain `Individual_ID` and `genoSample_ID` to match genotype plus covariates to correct for in the phenotype association. Column `Dx` (0 control 1 case) is optional, if present is used to build the reference set when computing T-scores. *Note: samples in genotype and phenotype matrix are matched based on covariate matrix*
- **Reactome Pathway annotation** (*--reactome_file*): .gmt file can be downloaded from https://reactome.org/download-data/ (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **GO Pathway annotation** (*--GOterms_file*): .RData file, can be obtained using *Annotate_GOterm_run.R*, each pathway is a entry in the list with `GOID` `Term` `Ontology` `geneIds` elemnets (provided in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/))
- **Custom pathway**: .RData file, similar to GO structure, each pathway is a list entry with `name` and `geneIds` elements. Available for WikiPathways (2019) in [refData](https://gitlab.mpcdf.mpg.de/luciat/castom-igex/-/tree/master/refData/)

## Workflow
### Predict gene expression
From previously trained PriLer tissue-specific model (Module 1), predict gene expression based on genotype-only dataset

*--InfoFold* is the folder with gene-snp distance matrix ENSEMBL_gene_SNP_2e+5_chr<>_matrix.mtx, *--outTrain_fold* is the folder with overall results from PriLer

#### Usage
```sh
./Priler_predictGeneExp_run.R \
    --genoDat_file \
    --covDat_file \
    --cis_thres (default = 200000) \
    --outFold \
    --outTrain_fold  \
    --InfoFold
```
*NOTE: can be splitted for subset of samples (depends on covDat_file), genes are NOT filtered*

The output includes (saved in *--outFold*):
- predictedExpression.txt.gz 
 
Based on data dimension, the next scripts are divided in two parts. If sample size <= 10,000 follow "Small dataset" part, otherwise "Large dataset".

***
***

### Small dataset: T-scores and Pathway-scores computation

Tscore_PathScore_diff_run.R

### Small dataset: Association with phenotype of T-score and pathways

pheno_association_smallData_run.R/pheno_association_smallData_customPath_run.R

### Small dataset: Meta-analysis across multiple cohorts

pheno_association_metaAnalysis_run.R/pheno_association_customPath_metaAnalysis_run.R

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
> Rscript PathwayScores_splitGenes_run.R --ncores (default 10) --input_file  --covDat_file  --outFold --split_tot (default 100) --reactome_file --GOterms_file --skip_reactome (default F)
- *input_file*: common path to .RData object predicted Tscores ({i}.RData part excluded)

The output includes:
- Pathway_Reactome/GO_scores.RData: pathway scores for each pathways and samples
- Pathway_Reactome/GO_pvalues.RData: pathway pvalues (t.test) for each pathways and samples
- PathwaySummaryPvalues_cases_Reactome/GO.txt: Fisher’s combined probability pvalue consider only cases
- PathwaySummaryPvalues_cases_Reactome/GO.txt: Fisher’s combined probability pvalue consider only controls

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


