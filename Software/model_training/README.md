# Training prediction model with PriLer

PriLer (CASTom-iGEx Module 1) is a pipeline (R based) that create model for predicting gene expression from genotype data integrating prior features on variants whose relevance is automatically learned in a machine learning set up. 

## Input Files
- **Gene expression matrix** (*--geneExp_file*): preprocessed gene expression (genes x samples). First column refers to gene names (ensembl annotation or HUGO nomenclature)
- **Genotype matrix** (*--genoDat_file*): dosages for each chromosome (compressed txt) without variants name/position (variants x samples). *NOTE: the file must end with chr<>_matrix.txt.gz*
- **Genotype info matrix** (*--VarInfo_file*): contains variants position, name and other info, must contain columns `CHR` and `POS` and match with Genotype matrix). *NOTE: the file must end with  chr<>.txt*
- **Covariate matrix** (*--covDat_file*): columns must contain `Individual_ID`, `genoSample_ID` and `RNASample_ID` to match genotype and gene expression plus covariates to be used in the regression model. Column `Dx` (0 control 1 case) is optional as well as it's usage. *Note: Samples in genotype and gene expression matrix are filtered based on covariate matrix*
- **Prior matrix** (*--priorDat_file*): prior information for variants (variants x prior features). It doesn't include variant name and MUST match genotype matrix. The columns can be binary (one-hot encoding for intersection) or continuous. 
- **List heritable genes** (*--geneList_file*): usually obtained from TWAS heritable analysis: list of genes, match external_gene_name or ensembl_gene_id. Set of heritable genes being regulated by cis-variants.
- **Gene annotation files of TSS and position** (*--biomartTSS_file --biomartGenePos_file*) obtained using *PrepareData_biomart_TSS.R* script. Possibility of recomputing or use provided fixed version

## Workflow
### Pre-processing:
Prepare files needed for the regression model analysis: 
annotate genes file using bioMart (possibility of recomputing or use the fixed version), compute snp-gene distance sparse matrix. If list of heritable genes not provided, all genes are annotated as not heritable.

*--outFold*, *--outFold_geneExp* and *--outFold_snps* are the path where general output, filtered gene expression and SNP annotation shoudl be saved respectively. 

*NOTE: consider only chromosomes 1-22*
```sh
./preProcessing_data_run.R \
	--geneExp_file \
	--geneList_file (default NULL) \
	--VarInfo_file \
	--cis_thres (default 200000) \
	--biomartGenePos_file (default NULL) \
	--biomartTSS_file (default NULL) \
	--outFold \
	--outFold_geneExp \
	--outFold_snps (default NULL)
```

The output includes:
-   RNAseq_filt.txt: gene expression + annotation table (genes x samples), saved in *--outFold_geneExp*
-   hg19_ENSEMBL_TSS_chr<>_matched.txt: gene annotation include column with heritable info, saved in *--outFold*
-   hg19_SNPs_chr<>_matched.txt: snp annotation, containts position and ID, saved in *--outFold_snps*
-   ENSEMBL_gene_SNP_2e+5_chr<>_matrix.mtx: sparse matrix (variants x genes) indicates the distance from the gene in bp with a default threshold of 200 kb, saved in *--outFold*

### Step 1:
Considering only heritable genes, compute elastic-net regression in a nested cross validation setting without prior information. The aim is to find the optimal alpha-lambda couple parameter for the outer loop. In addition, regression without prior is evaluated.

*--InfoFold* equivalent to *--outFold* in preProcessing and is the folder including ENSEMBL_gene_SNP_2e+5_chr<>_matrix.mtx data. *--curChrom* indicates considered choromosome in e.g. chr1

*NOTE: script is specific for a chromosome and can be used if no genes are heritable such that e-net is computed for all genes*
```sh
./PriLer_part1_run.R \
	--curChrom \
	--covDat_file \
	--genoDat_file \
	--geneExp_file \
	--outFold \
	--InfoFold \
	--functR ./Priler_functions.R \
	--ncores (default 20) \
	--seed_out (default 1234) \
	--seed_in (default 42) \
	--nfolds_in (default 5) \
	--nfolds_out (default 5) \
	--cis_thres (default 200000) \
	--Dx (default F)
```
The ouput includes (saved in *--outFold*):
-   optim_lambda_chr<>.txt/optim_alpha_chr<>.txt: optimal lambda/alpha parameter for each outer fold and gene (genes x outer folds)
-   resNoPrior_NestedCV_HeritableGenes_chr<>.RData/resNoPrior_NestedCV_AllGenes_chr<>.RData: R object containing    
	-	geneAnn: gene annotation
	-  	train: evaluation on train set for each outer folder
	-   test: evaluation on test set for each outer folder
	-   cor_comb_test: correlation original (cov.adjusted) vs predicted expression combing outer test folders 
	-	cor_comb_noadj_test: correlation original (not adjusted) vs predicted expression combing outer test folders
	-   beta_snps: regression coefficient for variants in each outer folder
	-   beta_cov: regression coefficient for covariates in each outer folders
	-   seed: seed to generate inner and outer partitions

### Step 2:
Considering only heritable genes, compute elastic-net regression in a nested cross-validation setting using prior information in order to find optimal E (scale for prior weights) parameter, alpha and lambda are obtained from the previous step. 

*--part1Res_fold* is the folder containing the output of step1, *--priorInf* is a vector of indeces indicating the columns in priorDat_file to be included.

*NOTE: The script is parallelized over given possible values of E parameter. Possible values for E cannot be chosen a prior but depends on the data.*

```sh
./PriLer_part2_run.R \
	--covDat_file \
	--genoDat_file \
	--geneExp_file \
	--InfoFold \
	--part1Res_fold \
	--priorDat_file \
	--priorInf (default 0)  \
	--ncores (default 10) \
	--functR ./Priler_functions.R \
	--cis_thres (default 200000) \
	--Dx (default F) \
	--maxIter (default 20) \
	--dThres  (default 0.001) \
	--convert_par (default 0.25) \
	--E_set \
	--outFold
```

The output includes (saved in *--outFold*):
 -   resE_allchr.RData:  R object with info of E parameter search for each E parameter, folder and interation (not further use, only kept to check/specific plots).
-   cv_train_Epar_allchr.txt/cv_test_Epar_allchr.txt: n.folds x n.Epar tested, sum of mean squared error (MSE) on train/test set for all the genes in the final iteration.
-   obj_Epar_cvtrain_allchr.RData/obj_Epar_cvtest_allchr.RData:  list for each E parameter and each folder of objective function on train/test sets (all iterations). 
-   evalf_Epar_cvtrain_allchr.RData/evalf_Epar_cvtest_allchr.RData  list for each E parameter and each folder of evaluation function (sum of MSE) on train/test set (all iterations). 
-   resPrior_EOpt(orEFixed)_NestedCV_HeritableGenes_allchr.RData  R object containing results for the optimal E parameter and the corresponding latest iteration. If minimum is reached as the latest E value (i.e. not reached in the interval given), then optimal E is chosen as the paramter that lead to objective function convergence in decreasing: |obj(step_i) - obj(step_i+1)|< 0.5
	-   geneAnn: gene annotation
	-   train_opt: evaluation on train set for each outer folder
	-   test_opt: evaluation on test set for each outer folder
	-   cor_comb_test_opt: correlation original vs predicted expression combing outer test folder
	-	cor_comb_test_noadj_opt: correlation original (not adjusted) vs predicted expression combining outer test folder    
	-   beta_snps_opt: regression coefficient for SNPs in each outer folder (divided by chr)
	-   beta_cov_opt: regression coefficient for covariates in each outer folder (divided by chr)
	-   weights_opt: n.prior features x n.outer folds, weights associated to each prior feature
	-	E_opt: value of optimal E paramter

#### E parameter calibration:
Start with 10 values (0.5, 1, 1.5, 2-8) and run part2. To evaluate if the range include the optimal value check average MSE in cv_test_Epar_allchr.txt. The first (few) E parameter(s) should give bigger errore with wide sd compared to the rest, afterward the error should decrease until convergence. Before convergence, there should be a value of E for which the average MSE is minimal.
If there is no "big gap", move the window to smaller value of E. If there is but no minimum for CV test that can be found, add additional intermediate points. In case a minimum cannot be found, the algorithm chose the first value of convergence, although not optimal.

### Step 3:
Considering only heritable genes, first find optimal alpha and lambda parameter on the entire set (single cross validation) and evaluate total results without prior. Second, use alpha-lambda pairs found and the optimal E parameter (step 2) in the elastic-net with prior information setting and evaluate the results.

*--part2Res_fold* is the folder containing the output of step2.

```sh
./PriLer_part3_run.R \
	--covDat_file \
	--genoDat_file \
	--geneExp_file \
	--InfoFold \
	--part2Res_fold \
	--priorDat_file \
	--priorInf (default 0) \
	--ncores (default 10) \
	--functR ./Priler_functions.R \
	--cis_thres (default 200000) \
	--Dx (default F) \
	--maxIter (default 20) \
	--dThres (default 0.001) \
	--convert_par (default 0.25) \
	--seed (default 4321) \
	--nfolds (default 5) \
	--outFold
```

The output includes (saved in *--outFold*):
-   resNoPrior_HeritableGenes_allchr.RData: results without prior information:
	-   geneAnn: gene annotation
	-   tot: evaluation on entire set
	-   beta_snps: regression coefficient for SNPs, divided by chr
	-   beta_cov: regression coefficient for covariates, divided by chr
	-   seed: seed to generate folds
-   resPrior_EOpt(orEFixed)_Iteration_HeritableGenes_allchr.RData: R iteration results for prior convergence
	-   Eopt: optimal (or convergence) E parameter
	-   pWeight: n.iterations x n.prior features, prior weights at each iteration
	-   errComp: n.iterations x 3, sum MSE - penalty for beta - penalty for weights at each iteration    
	-   nCount: n.iterations x n. genes, number of variants regulating a gene
	-   dev/dev_geno/dev_cov/dev_genocov: n.iterations x n. genes, total deviance (R2), deviance explained by only genotype, deviance explained by only covariates, deviance explained by the interaction covariance/genotype
	-	cor/cor_pval:  n.iterations x n. genes, correlation original (cov.adjusted) vs predicted expression and corresponding p-value
	-	cor_noadj/cor_noadj_pval:  n.iterations x n. genes, correlation original (not adjusted) vs predicted expression and corresponding p-value
	-   MSE: n.iterations x n. genes, mean squared error 
	-   beta: for each chr n.variables x n. genes, regression coefficient for both variants and covariates at each step
	- 	obj: objective function at each step
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

*--part1Res_fold --part2Res_fold --part3Res_fold* are the folders containing output of step1, step2 and step3.

*NOTE: script is specific for a chromosome, seed and n. folds for nested CV and single CV are the same as used in the previous steps*

```sh
./PriLer_part4_run.R \
	--curChrom \
	--covDat_file \
	--genoDat_file \
	--geneExp_file \
	--ncores (default 20) \
	--outFold \
	--InfoFold \
	--functR ./PriLer_functions_run.R \
	--priorDat_file \
	--priorInf (default 0) \
	--part1Res_fold \
	--part2Res_fold \
	--part3Res_fold \
	--cis_thres (default 200000) \
	--Dx (default F) \
	--convert_par (default 0.25)
```

The output includes (saved in *--outFold*):
-   resNoPrior_NestedCV_NotHeritableGenes_chr<>.RData: R object as in PriLer_part1_run.R
-   resPrior_NestedCV_NotHeritableGenes_chr<>.RData: R object as above, uses prior in the regression
-   resNoPrior_NotHeritableGenes_chr<>.RData: R object as in PriLer_part3_run.R
-   resPrior_NotHeritableGenes_chr<>.RData: R object as above, uses prior in the regression

### Combine Output
The last step combines regression results from nested CV and total set  regression. In addition, results are separated for with and without setting to allow for a straightforward comparison. The regression coefficients are saved in a .RData object, divided by chromosomes. Finally, plot of the entire pipeline are produced. 

*--part1Res_fold --part2Res_fold --part3Res_fold --part4Res_fold* are the folders containing output of step1, step2, step3 and step4.
 
```sh
./PriLer_finalOutput_run.R \
	--covDat_file \
	--outFold \
	--InfoFold \
	--functR ./PriLer_functions_run.R \ 
	--priorDat_file \
	--priorInf (default 0) \
	--part1Res_fold \
	--part2Res_fold \
	--part3Res_fold \
	--part4Res_fold \
	--cis_thres (default 200000) \
	--Dx (default F) \
	--convert_par (default 0.25)
```

The output includes (saved in *--outFold*):
-   resNoPrior_regEval_allchr.txt and resPrior_regEval_allchr.txt regression evaluation for each gene + gene annotation
-   resNoPrior_regCoeffCov_allchr.txt and resPrior_regCoeffCov_allchr.txt regression coefficients for covariates
-   resNoPrior_regCoeffSnps_allchr.RData and resPrior_regCoeffSnps_allchr.RData regression coefficient for variants divided by chr (sparse matrix)
-   plots for model evaluation and quality checks
'NoPrior' refers to elastic-net, 'Prior' refers to PriLer results.

### Note on quality check
- *plot_devtest_compareWithPrior_increase.png* should indicate a percentage of increase > 50%
- *plot_CVmse_TEST_Eopt.png* should decrease at each iteration (at least for some folds)
- *plot_numberSNPs_elVSelP_gene.png* n. of reg-SNPs should be lower in PriLer with a decrease of "without prior" component
- *plot_iteration_weigths_priorFeat.png* when relevant, prior feature weigths iteratively increase 

