# Training prediction model with PriLer

PriLer is a command-line tool that create model for predicting gene expression from genotype data using prior features on variants whose relevance is automatically learned in a machine learning set up. 

## Requirements
To run PriLer the following R packages are required:

- argparse 
- Matrix
- biomaRt
- glmnet
- parallel
- doParallel
- nloptr
- bigmemory
- ggplot2
- gridExtra
- RColorBrewer
- ggsci
- ggExtra

## Input Files
- **Gene expression matrix**:  preprocessed gene expression. Genes on the rows, samples on the columns. First column refers to gene names (ensembl annotation or HUGO nomenclature)
- **Genotype matrix**: dosages for each chromosome (compressed txt) without variants name/position. Variants on the rows, samples on the columns.  *NOTE: the file must end with chr<>_matrix.txt.gz*
- **Genotype info matrix**: contains variants position, name and other info, must contain columns CHR and POS and match with Genotype matrix). *NOTE: the file must end with  chr<>.txt*
- **Covariate matrix**: columns must contain Individual_ID, genoSample_ID and RNASample_ID to match genotype and gene expression plus covariates to be used in the regression model. Column Dx (0 control 1 case) is optional as well as it's usage. *Note: samples in genotype and gene expression matrix are filtered based on covariate matrix*
- **Prior matrix**: prior information for variants. Variants on the rows, features on the columns. It doesn't include variant name and MUST match genotype matrix. The columns can be binary (only 2 values 0 no prior and 1 prior) or continuous. 
- **List heritable genes**: usually obtained from TWAS heritable analysis: list of genes, match external_gene_name or ensembl_gene_id. Reliable set of genes being regulated by cis-variants.
- **Gene annotation files of TSS and position** obtained using *PrepareData_biomart_TSS.R* script. Possibility of recomputing or use a fixed version

## Workflow
### Pre-processing:
Prepare files needed for the regression model analysis: 
annotate genes file using bioMart (possibility of recomputing or use the fixed version), compute snp-gene distance sparse matrix. *NOTE consider only chromosomes 1-22*
#### Usage
>Rscript preProcessing_data_run.R --geneExp_file --geneList_file --VarInfo_file --cis_thres (default 200000) --biomartGenePos_file (default NA) --biomartTSS_file (default NA) --outFold --outFold_snps (default NA)

The output includes:
-   RNAseq_filt.txt:  n.genes x n.samples gene expression + annotation table
-   hg19_ENSEMBL_TSS_chr<>_matched.txt: gene annotation include column with heritable info
-   hg19_SNPs_chr<>_matched.txt: snp annotation, containts position and ID
-   ENSEMBL_gene_SNP_2e+5_chr<>_matrix.mtx: sparse matrix n.variants x n.genes, indicates the distance from the gene in bp with a default threshold of 200 kb 

### Step 1:
Considering only heritable genes, compute elastic-net regression in a nested cross validation setting without prior information. The aim is to find the optimal alpha-lambda couple parameter for the outer loop. In addition, regression without prior is evaluated. 
*NOTE: script is specific fro a chromosome*
#### Usage
>Rscript ElNet_withPrior_part1_run.R --curChrom --covDat_file --genoDat_file --geneExp_file --ncores (default 20) --outFold --InfoFold --functR ElNet_withPrior_functions_run.R --seed_out (default 1234) --seed_in (default 42) --nfolds_in (default 5) --nfolds_out (default 5) --cis_thres (default 200000) --Dx (default F)

The ouput includes:
-   optim_lambda_chr<>.txt/optim_alpha_chr<>.txt: n.genes x n. outer folds optimal lambda/alpha parameter for each outer fold and gene
-   resNoPrior_NestedCV_HeritableGenes_chr<>.RData: R object containing    
	- geneAnn: gene annotation
	-   train: evaluation on train set for each outer folder
	-   test: evaluation on test set for each outer folder
	-   cor_comb_test: correlation true vs predicted expression combing outer test folders    
	-   beta_snps: regression coefficient for variants in each outer folder
	-   beta_cov: regression coefficient for covariates in each outer folders
	-   seed: seed to generate inner and outer partitions

### Step 2:
Considering only heritable genes, compute elastic-net regression in a nested cross-validation setting using prior information in order to find optimal E (scale for prior weights) parameter, alpha and lambda are obtained from the previous step. *NOTE: Possible values for E cannot be chosen a prior but depends on the data*.
#### Usage
> Rscript ElNet_withPrior_part2_run.R --covDat_file --genoDat_file --geneExp_file --InfoFold --part1Res_fold --priorDat_file --priorInf (default 0)  --ncores (default 10) --functR ElNet_withPrior_functions_run.R --cis_thres (default 200000)  --Dx (default F)  --maxIter 20 --dThres  0.001  --convert_par (default 0.25) --E_set --outFold

The output includes:
 -   resE_allchr.RData:  R object with info of E parameter search for each E parameter, folder and interation. (not further use, only kept to check/specific plots)
-   cv_train_Epar_allchr.txt/cv_test_Epar_allchr.txt: n.folds x n. Epar tested, sum of mean squared error (MSE) on train/test set for all the genes in the final iteration
-   obj_Epar_cvtrain_allchr.RData/obj_Epar_cvtest_allchr.RData:  list for each E parameter and each folder of objective function on train/test sets (all iterations) 
-   evalf_Epar_cvtrain_allchr.RData/evalf_Epar_cvtest_allchr.RData  list for each E parameter and each folder of evaluation function (sum of MSE) on train/test set (all iterations) 
-   resPrior_EOpt_NestedCV_HeritableGenes_allchr.RData  R object containing results for the optimal E parameter and the corresponding latest iteration
	-   geneAnn: gene annotation
	-   train_opt: evaluation on train set for each outer folder
	-   test_opt: evaluation on test set for each outer folder
	-   cor_comb_test_opt: correlation true vs predicted expression combing outer test folder    
	-   beta_snps_opt: regression coefficient for SNPs in each outer folder (divided by chr)
	-   beta_cov_opt: regression coefficient for covariates in each outer folder (divided by chr)
	-   weights_opt: n.prior features x n.outer folds, weights associated to each prior feature

### Step 3:
Considering only heritable genes, first find optimal alpha and lambda parameter on the entire set (single cross validation) and evaluate total results without prior. Second, use alpha-lambda pairs found and the optimal E parameter in the elastic-net with prior information setting and evaluate the results. 

#### Usage
> Rscript ElNet_withPrior_part3_run.R --covDat_file --genoDat_file --geneExp_file --InfoFold --part2Res_fold --priorDat_file --priorInf (default 0)  --ncores (default 10) --functR ElNet_withPrior_functions_run.R  --cis_thres (default 200000)  --Dx (default F)  --maxIter 20 --dThres  0.001  --convert_par (default 0.25) --outFold --seed (default 4321) --nfolds (default 5)

The output includes:
-   resNoPrior_HeritableGenes_allchr.RData: results without prior information:
	-   geneAnn: gene annotation
	-   tot: evaluation on entire set
	-   beta_snps: regression coefficient for SNPs, divided by chr
	-   beta_cov: regression coefficient for covariates, divided by chr
	-   seed: seed to generate folds
-   resPrior_EOpt(orEFixed)_Iteration_HeritableGenes_allchr.RData: R iteration results for prior convergence
	-   Eopt: optimal (or convergence) E parameter
	-   pWeight: n.iterations x n.prior features, prior weights at each iteration
	-   errComp: n.Iterations x 3, sum MSE - penalty for beta - penalty for weights at each iteration    
	-   nCount: n.Iterations x n. genes, number of variants regulating a gene
	-   dev/dev_geno/dev_cov/dev_genocov: n.Iterations x n. genes, total deviance (R2), deviance explained by only genotype, deviance explained by only covariates, deviance explained by the interaction covariance/genotype
	- cor/cor_pval:  n.Iterations x n. genes, correlation true vs predicted expression and corresponding p-value
	-   MSE: n.Iterations x n. genes, mean squared error 
	-   beta: for each chr n.variables x n. genes, regression coefficient for both variants and covariates at each step
	-  obj: objective function at each step
	- 
-   resPrior_EOpt(orEFixed)_HeritableGenes_allchr.RData: results with prior information (final iteration)	    
	-   geneAnn: gene annotation
	-   tot: evaluation on entire set
	-   beta_snps: regression coefficient for SNPs, divided by chr
	-   beta_cov: regression coefficient for covariates, divided by chr
	-   Eopt: optimal (or convergence) E parameter
	-   weights: final weights for prior features
	-   prior_coeff: for each chromosome, prior coefficient for each variants

### Step 4:
Considering only not heritable genes, fix prior coefficients found in the previous steps (step 4 in the total set and step 2 in the outer folds for nested CV). First, compute elastic-net regression in a nested CV setting without prior information, find optimal alpha-lambda pairs and evaluate the regression without prior. Second, use optimal alpha-lambda and prior coefficients to evaluate elastic-net regression in a nested CV setting with prior information. Third, find optimal alpha-lambda parameters on the entire set (single CV), evaluate total results with and without prior.
*NOTE: script is specific for a chromosome, seed and n. folds for nested CV and single CV are the same as used in the previous steps*

#### Usage
> Rscript ElNet_withPrior_part4_run.R --curChrom --covDat_file --genoDat_file --geneExp_file --ncores (default 20) --outFold --InfoFold --functR ElNet_withPrior_functions_run.R --priorDat_file --priorInf (default 0) --part1Res_fold --part2Res_fold --part3Res_fold --cis_thres (default 200000) --Dx (default F) --convert_par (default 0.25)

The output includes:
-   resNoPrior_NestedCV_NotHeritableGenes_chr<>.RData: R object as in ElNet_withPrior_part1_run.R
-   resPrior_NestedCV_NotHeritableGenes_chr<>.RData: R object as above, uses prior in the regression
-   resNoPrior_NotHeritableGenes_chr<>.RData: R object as in ElNet_withPrior_part3_run.R
-   resPrior_NotHeritableGenes_chr<>.RData: R object as above, uses prior in the regression

### Combine Output
The last step combines regression results from nested CV and total set  regression. In addition, results are separated for with and without setting to allow for a straightforward comparison. The regression coefficients are saved in a .RData object, divided by chromosomes. Finally, plot of the entire pipeline are produced. 
#### Usage 
>Rscript ElNet_withPrior_finalOutput_run.R --covDat_file --outFold --InfoFold --functR ElNet_withPrior_functions_run.R --priorDat_file --priorInf (default 0) --part1Res_fold --part2Res_fold --part3Res_fold --part4Res_fold --cis_thres (default 200000) --Dx (default F) --convert_par (default 0.25)

The output includes:
-   resNoPrior_regEval_allchr.txt and resPrior_regEval_allchr.txt regression evaluation for each gene + gene annotation
-   resNoPrior_regCoeffCov_allchr.txt and resPrior_regCoeffCov_allchr.txt regression coefficients for covariates
-   resNoPrior_regCoeffSnps_allchr.RData and resPrior_regCoeffSnps_allchr.RData regression coefficient for variants divided by chr (sparse matrix)
-   plots
