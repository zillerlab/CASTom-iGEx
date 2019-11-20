#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(parallel))

parser <- ArgumentParser(description="part4: regression for not heritable genes, with and without found prior, divided by chr")

parser$add_argument("--curChrom", type = "character", help = "chromosome")
parser$add_argument("--genoDat_file", type = "character", help = "genotype data, common path (ending chr specific)")
parser$add_argument("--geneExp_file", type = "character", help = "RNA expression, complete path")
parser$add_argument("--covDat_file", type = "character", help = "covariance file, complete path")
parser$add_argument("--InfoFold", type = "character", help = "path to fold with additional info (e.g. snp-gene dist matrix)")
parser$add_argument("--priorDat_file", type = "character", help = "prior matrix file, common path (ending chr specific)")
parser$add_argument("--priorInf", type="integer", nargs = '*' , default = 0, help = "index prior feature to be used, if 0 all column are used")
parser$add_argument("--ncores", type="integer", default = 20, help = "n.of cores for the parallelization")
parser$add_argument("--functR", type="character", help = "Rscript with functions to be used, complete path")
parser$add_argument("--part1Res_fold", type = "character", help = "path folder reuslt part1, seed and nfolds for nested CV needed")
parser$add_argument("--part2Res_fold", type = "character", help = "path folder reuslt part2, prior for each outer folder needed")
parser$add_argument("--part3Res_fold", type = "character", help = "path folder reuslt part3, seed and nfolds for CV needed, total prior")
parser$add_argument("--cis_thres", type = "integer", default = 200000, help = "window (in bp) to compute distance from gene and snps")
parser$add_argument("--Dx", type = "logical", default = FALSE, help = "if true, Dx included as covariate")
parser$add_argument("--convert_par", type="double",  default = 0.25, help = "")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
curChrom <- args$curChrom
genoDat_file <- args$genoDat_file
geneExp_file <- args$geneExp_file
covDat_file <- args$covDat_file
InfoFold <- args$InfoFold
priorDat_file <- args$priorDat_file
priorInf <- args$priorInf
ncores <- args$ncores
functR <- args$functR
part1Res_fold <- args$part1Res_fold
part2Res_fold <- args$part2Res_fold
part3Res_fold <- args$part3Res_fold
cis_thres <- args$cis_thres
Dx <- args$Dx
convert_par <- args$convert_par
outFold <- args$outFold


#################################################################
# curChrom <- 'chr21'
# covDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/covariateMatrix_Control50.txt'
# genoDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/Genotyping_data/Genotype_dosage_'
# geneExp_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/RNAseq_data/EXCLUDE_ANCESTRY_SVA/RNAseq_filt.txt'
# ncores <- 10
# cis_thres = 200000
# outFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/'
# InfoFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/'
# functR <- '/ziller/lucia/eQTL_PROJECT_CMC/RSCRIPTS/SCRIPTS_v1/ElNet_withPrior_functions_run.R'
# Dx = F
# part1Res_fold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/200kb/'
# part2Res_fold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/200kb/'
# part3Res_fold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/200kb/'
#################################################################

source(functR)
print(paste('Include Dx in covariates:',Dx))


#####################
### load datasets ###
#####################
## covariates matrix, contains sample annotation
covDat <- read.table(covDat_file, header = T, sep = '\t', stringsAsFactors = F, check.names=F)
sampleAnn <- covDat[, colnames(covDat) %in% c('Individual_ID', 'genoSample_ID', 'RNASample_ID')]

if(!Dx){col_exclude <- c('Individual_ID', 'genoSample_ID', 'RNASample_ID', 'Dx')}
covDat <- covDat[, !colnames(covDat) %in% col_exclude]

## RNA-seq
expDat <- read.table(geneExp_file, header = T, sep = '\t', stringsAsFactors = F, check.names=F)
gene_ann <- expDat[expDat$chrom == curChrom, colnames(expDat) %in% c('type',	'chrom', 'TSS_start','TSS_end','name','start_position','end_position','ensembl_gene_id','external_gene_name')]
expDat <- expDat[expDat$chrom == curChrom, !colnames(expDat) %in% c('type',	'chrom', 'TSS_start','TSS_end','name','start_position','end_position','ensembl_gene_id','external_gene_name')]
# order based on sample ann
id <- unname(sapply(sampleAnn$RNASample_ID, function(x) which(x == colnames(expDat))))
expDat <- expDat[, id]
# consider only heritable genes and save the id
id_nother <- which(gene_ann$type == 'not_heritable')
expDat <- expDat[id_nother, ]
gene_ann <- gene_ann[id_nother,]

## genotype
genDat <- read.table(gzfile(sprintf('%s%s_matrix.txt.gz', genoDat_file, curChrom)), header = T, sep = '\t', check.names=F)
# order based on sample ann
id <- unname(sapply(sampleAnn$genoSample_ID, function(x) which(x == colnames(genDat))))
genDat <- genDat[, id]

## gene-snp distance matrix
cis_ann <- paste0(cis_thres/10^5, 'e+05')
print(cis_ann)
geneSnpDist <- readMM(sprintf('%sENSEMBL_gene_SNP_%s_%s_matrix.mtx', InfoFold, cis_ann, curChrom))
# only heritable genes
geneSnpDist <- geneSnpDist[, id_nother]
covIndex <- c((nrow(geneSnpDist)+1):(nrow(geneSnpDist)+ncol(covDat)))

pDat <- read.table(gzfile(sprintf('%s%s.txt.gz', priorDat_file, curChrom)), header = T, sep = '\t')
pNames <- colnames(pDat)[priorInf]
pDat <- as.matrix(pDat[,priorInf]) # only cell type of interest
# pNames <- colnames(pDat)
print(pNames)
priorMat <- pDat*convert_par
print(head(priorMat))

# # load prior for all chromosomes, used to normalize binary values
# all_Chroms <- paste0('chr', 1:22)
# pDat <- list()
# for(i in 1:length(all_Chroms)){
#   
#   # genotype
#   chr <- all_Chroms[i]
#   print(chr)
#   
#   ## prior matrix
#   pDat[[i]] <- read.table(gzfile(sprintf('%s%s.txt.gz', priorDat_file, chr)), header = T, sep = '\t')
#   if(length(priorInf) == 1 & priorInf[1] == 0){priorInf <- 1:ncol(pDat[[i]])}
#   pDat[[i]] <- as.matrix(pDat[[i]][,priorInf]) # only cell type of interest
#   
# }
# 
# # normalize epi priors
# type_val <- apply(pDat[[1]], 2, unique)
# continous_pos <- !sapply(type_val, function(x) all(x %in% c(0,1)) & length(x)==2) 
# tmp <- do.call(rbind, pDat)
# pNames <- colnames(tmp)

# print(pNames)

#if(length(which(!continous_pos))>1){
#  N_k <- apply(tmp[,!continous_pos], 2, function(x) length(which(x != 0)))
#  pDat <- lapply(pDat, function(x) cbind(x[,continous_pos], sweep(x[,!continous_pos], 2, N_k/min(N_k), '/')))
#  # pDat[, !gwas_pos] <- sweep(pDat[,!gwas_pos], 2, N_k/min(N_k), '/')
#}

# N_k <- apply(tmp, 2, function(x) mean(x))
# pDat <- lapply(pDat, function(x)  sweep(x, 2, N_k/min(N_k), '/'))


chr <- as.numeric(strsplit(curChrom, split = 'chr')[[1]][2])
# priorMat <- pDat[[chr]]


#############################
### load previous results ###
#############################

# seed and nfolds for nested CV
resNoPrior_Her_nCV <- get(load(sprintf('%s/resNoPrior_NestedCV_HeritableGenes_%s.RData', part1Res_fold, curChrom)))
seed_in <- resNoPrior_Her_nCV$seed$value[resNoPrior_Her_nCV$seed$type == 'inner']
nfolds_in <- resNoPrior_Her_nCV$seed$nfold[resNoPrior_Her_nCV$seed$type == 'inner']
seed_out <- resNoPrior_Her_nCV$seed$value[resNoPrior_Her_nCV$seed$type == 'outer']
nfolds_out <- resNoPrior_Her_nCV$seed$nfold[resNoPrior_Her_nCV$seed$type == 'outer']

# prior for nested CV
resPrior_Her_nCV <- get(load(sprintf('%s/resPrior_EOpt_NestedCV_HeritableGenes_allchr.RData',part2Res_fold)))
priorCoeff_nCV <- apply(resPrior_Her_nCV$weights_opt, 2, function(x) getPrior(priorMat, x))

# total prior and seed/nfold CV
err_test <- read.table(sprintf('%s/cv_test_Epar_allchr.txt', part2Res_fold), header = T)
Epar <- sapply(colnames(err_test), function(x) as.numeric(strsplit(x, 'E_h.')[[1]][2]))
names(Epar) <- NULL
id_con <- which(abs(diff(colMeans(err_test))) < 0.5)[1] + 1
id_opt <- which.min(colMeans(err_test)) 
names(id_opt) <- NULL
names(id_con) <- NULL
if(id_opt != length(Epar)){
  id_con <- id_opt
}

if(id_con == id_opt){
  resPrior_Her <- get(load(sprintf('%s/resPrior_EOpt_HeritableGenes_allchr.RData', part3Res_fold)))
}else{
  resPrior_Her <- get(load(sprintf('%s/resPrior_EFixed%.2f_HeritableGenes_allchr.RData', part3Res_fold, Epar[id_con])))
}
priorCoeff <- resPrior_Her$prior_coeff[[chr]]

resNoPrior_Her <- get(load(sprintf('%s/resNoPrior_HeritableGenes_allchr.RData', part3Res_fold)))
seed <- resNoPrior_Her$seed$value
nfolds <- resNoPrior_Her$seed$nfold

#####################################################################################################################################
# input values
alpha_h <- seq(0.1, 0.9, by = 0.1) # exclude lasso because it choses only 1 varaible different from zero between the correlated ones 
# no need to standardize!

genDat <- t(genDat)
expDat <- t(expDat)

##
N <- ncol(expDat) # n. of genes
M <- nrow(expDat) # n. of samples
P <- ncol(genDat) # n. of SNPs
print(paste('genes =', N, ' variants =', P, ' samples =', M))
##

###########################################
#### nested CV with and without prior #####
###########################################

# divide in fold 
# use gradient descent algorithm
set.seed(seed_out)
Folds <- generateCVRuns(nfold = nfolds_out, labels = 1:M, ntimes = 1)[[1]]

# use iterative coordinate descent
clust <- makeCluster(ncores)
clusterExport(clust, "geneSnpDist")
clusterExport(clust, "covIndex")
clusterExport(clust, "covDat")
clusterExport(clust, "expDat")
clusterExport(clust, "genDat")
clusterEvalQ(clust,library(glmnet))
clusterExport(clust,'generateCVRuns')

lambdaVec <- vector(mode = 'list', length = nfolds_out)
alphaVec <- vector(mode = 'list', length = nfolds_out)

res_nop <- vector(mode = 'list', length = nfolds_out)
test_res_nop <- vector(mode = 'list', length = nfolds_out)
betaSnpsMatrix_nop <- vector(mode = 'list', length = nfolds_out)
betaCovMatrix_nop <- vector(mode = 'list', length = nfolds_out)

res_p <- vector(mode = 'list', length = nfolds_out)
test_res_p <- vector(mode = 'list', length = nfolds_out)
betaSnpsMatrix_p <- vector(mode = 'list', length = nfolds_out)
betaCovMatrix_p <- vector(mode = 'list', length = nfolds_out)


# not prior info
for(k in 1:nfolds_out){
  
  print(k)
  
  #### whitout prior ####
  system.time(res_nop[[k]] <- parLapply(clust, seq(1,N), expPrediction_cv_noPrior_fold, alpha = alpha_h, 
                                    seed = seed_in, nfolds = nfolds_in, fold = setdiff(1:M, Folds[[k]])))
  res_nop[[k]] <- do.call(rbind, res_nop[[k]])
  colnames(res_nop[[k]]) <- c('gene', 'idBeta', 'Beta', 'lambda', 'alpha', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')
  res_nop[[k]] <- as.data.frame(res_nop[[k]])
  print(head(res_nop[[k]]))
  
  lambdaVec[[k]] <- res_nop[[k]][!duplicated(res_nop[[k]][,1]),4]
  lambdaVec[[k]][is.na(lambdaVec[[k]])] <- 0
  alphaVec[[k]] <- res_nop[[k]][!duplicated(res_nop[[k]][,1]),5]
  alphaVec[[k]][is.na(alphaVec[[k]])] <- 0
  
  # save beta in sparse matrix format
  tmp <- res_nop[[k]][!is.na(res_nop[[k]][,2]),]
  modelMatrix <- sparseMatrix(i=tmp[,2],j=tmp[,1],x=tmp[,3],dims=c(P + ncol(covDat) + 1, N), giveCsparse=T)
  
  betaSnpsMatrix_nop[[k]] <- modelMatrix[1:P,]
  betaCovMatrix_nop[[k]] <- modelMatrix[(P+1):nrow(modelMatrix),]
  
  test_res_nop[[k]] <- deviance_funct(genDat, expDat, covDat, modelMatrix, Folds[[k]])
  colnames(test_res_nop[[k]]) <- c('dev', 'dev_geno','dev_cov', 'dev_geno_cov', 'cor', 'cor_pval')
  test_res_nop[[k]] <- as.data.frame(test_res_nop[[k]])
  
  #### with prior, fixed alpha and lambda ####
  res_p[[k]] <- parLapply(clust, seq(1,N), expPrediction_Prior_fold, prior = as.vector(priorCoeff_nCV[,k]), lambda =  lambdaVec[[k]], 
                                alpha = alphaVec[[k]], fold = setdiff(1:M, Folds[[k]]))
  
  res_p[[k]] <- do.call(rbind, res_p[[k]])
  colnames(res_p[[k]]) <- c('gene', 'idBeta', 'Beta', 'lambda', 'alpha', 'SquErr', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')
  res_p[[k]] <- as.data.frame(res_p[[k]])
  print(head(res_p[[k]]))
  
  # save beta in sparse matrix format
  tmp <- res_p[[k]][!is.na(res_p[[k]][,2]),]
  modelMatrix <- sparseMatrix(i=tmp[,2],j=tmp[,1],x=tmp[,3],dims=c(P + ncol(covDat) + 1, N), giveCsparse=T)
  
  betaSnpsMatrix_p[[k]] <- modelMatrix[1:P,]
  betaCovMatrix_p[[k]] <- modelMatrix[(P+1):nrow(modelMatrix),]
  
  test_res_p[[k]] <- deviance_funct(genDat, expDat, covDat, modelMatrix, Folds[[k]])
  colnames(test_res_p[[k]]) <- c('dev', 'dev_geno','dev_cov', 'dev_geno_cov', 'cor', 'cor_pval')
  test_res_p[[k]] <- as.data.frame(test_res_p[[k]])
  
  print(paste('regression computed fold', k))
  
  
}

predTest_geno_comb_nop <- lapply(1:nfolds_out, function(x) as.matrix(t(betaSnpsMatrix_nop[[k]]) %*% t(genDat[Folds[[x]],])))
predTest_geno_comb_nop <- do.call(cbind, predTest_geno_comb_nop)
predTest_cov_comb_nop <- lapply(1:nfolds_out, function(x) as.matrix(t(betaCovMatrix_nop[[k]]) %*% t(cbind(covDat[Folds[[x]],], rep(1,length(Folds[[x]]))))))

adjExpTest_comb_nop <- lapply(1:nfolds_out, function(x) t(expDat[Folds[[x]],]) - predTest_cov_comb_nop[[x]])
adjExpTest_comb_nop <- do.call(cbind, adjExpTest_comb_nop)

cor_comb_test_nop <- data.frame(cor = rep(0,N), cor_pval = rep(0,N))
for(g in 1:N){
  cor_comb_test_nop$cor[g] <- cor.test(adjExpTest_comb_nop[g,],predTest_geno_comb_nop[g,])$estimate
  cor_comb_test_nop$cor_pval[g] <- cor.test(adjExpTest_comb_nop[g,],predTest_geno_comb_nop[g,])$p.value
}


predTest_geno_comb_p <- lapply(1:nfolds_out, function(x) as.matrix(t(betaSnpsMatrix_p[[k]]) %*% t(genDat[Folds[[x]],])))
predTest_geno_comb_p <- do.call(cbind, predTest_geno_comb_p)
predTest_cov_comb_p <- lapply(1:nfolds_out, function(x) as.matrix(t(betaCovMatrix_p[[k]]) %*% t(cbind(covDat[Folds[[x]],], rep(1,length(Folds[[x]]))))))

adjExpTest_comb_p <- lapply(1:nfolds_out, function(x) t(expDat[Folds[[x]],]) - predTest_cov_comb_p[[x]])
adjExpTest_comb_p <- do.call(cbind, adjExpTest_comb_p)

cor_comb_test_p <- data.frame(cor = rep(0,N), cor_pval = rep(0,N))
for(g in 1:N){
  cor_comb_test_p$cor[g] <- cor.test(adjExpTest_comb_p[g,],predTest_geno_comb_p[g,])$estimate
  cor_comb_test_p$cor_pval[g] <- cor.test(adjExpTest_comb_p[g,],predTest_geno_comb_p[g,])$p.value
}


train_res_nop <- lapply(res_nop, function(x) x[!duplicated(x$gene), colnames(x) %in% c('lambda','alpha', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')])
train_res_p <- lapply(res_p, function(x) x[!duplicated(x$gene), colnames(x) %in% c('dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')])

res_tot_nop <- list(geneAnn = gene_ann, train = train_res_nop, test = test_res_nop, cor_comb_test = cor_comb_test_nop,
                beta_snps = betaSnpsMatrix_nop, beta_cov = betaCovMatrix_nop, 
                seed = data.frame(type = c('inner', 'outer'), value = c(seed_in, seed_out), nfold = c(nfolds_in, nfolds_out)))

res_tot_p <- list(geneAnn = gene_ann, train = train_res_p, test = test_res_p, cor_comb_test = cor_comb_test_p,
                    beta_snps = betaSnpsMatrix_p, beta_cov = betaCovMatrix_p, 
                    seed = data.frame(type = c('inner', 'outer'), value = c(seed_in, seed_out), nfold = c(nfolds_in, nfolds_out)))

# store results
save(res_tot_nop, file = sprintf('%s/resNoPrior_NestedCV_NotHeritableGenes_%s.RData', outFold, curChrom))
save(res_tot_p, file = sprintf('%s/resPrior_NestedCV_NotHeritableGenes_%s.RData', outFold,  curChrom))


#####################################################
##### res on all samples with and without prior #####
#####################################################

res_nop_alls <-  parLapply(clust, seq(1,N), expPrediction_cv_noPrior_fold, alpha = alpha_h, seed = seed, nfolds = nfolds, fold = 1:M)
res_nop_alls <- do.call(rbind, res_nop_alls)
colnames(res_nop_alls) <- c('gene', 'idBeta', 'Beta', 'lambda', 'alpha', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')
res_nop_alls <- as.data.frame(res_nop_alls)

tmp <- res_nop_alls[!is.na(res_nop_alls[,2]),]
modelMatrix <- sparseMatrix(i=tmp[,2],j=tmp[,1],x=tmp[,3],dims=c(P + ncol(covDat) + 1, N), giveCsparse=T)

betaSnpsMatrix_nop_tot <- modelMatrix[1:P,]
betaCovMatrix_nop_tot <- modelMatrix[(P+1):nrow(modelMatrix),]

tot_nop_alls <- res_nop_alls[!duplicated(res_nop_alls$gene), colnames(res_nop_alls) %in% c('lambda', 'alpha', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')]

res_tot_nop <- list(geneAnn = gene_ann, tot = tot_nop_alls, beta_snps = betaSnpsMatrix_nop_tot, beta_cov = betaCovMatrix_nop_tot, 
                    seed = data.frame(type ='tot', value = seed, nfold = nfolds)) 

save(res_tot_nop, file = sprintf('%s/resNoPrior_NotHeritableGenes_%s.RData', outFold,  curChrom))

## 
res_p_alls <- parLapply(clust, seq(1,N), expPrediction_Prior_fold, prior = priorCoeff, lambda =  tot_nop_alls$lambda, 
                        alpha = tot_nop_alls$alpha, fold = 1:M)

res_p_alls <- do.call(rbind, res_p_alls)
colnames(res_p_alls) <- c('gene', 'idBeta', 'Beta', 'lambda', 'alpha', 'SquErr', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')
res_p_alls <- as.data.frame(res_p_alls)

tmp <- res_p_alls[!is.na(res_p_alls[,2]),]
modelMatrix <- sparseMatrix(i=tmp[,2],j=tmp[,1],x=tmp[,3],dims=c(P + ncol(covDat) + 1, N), giveCsparse=T)

betaSnpsMatrix_p_tot <- modelMatrix[1:P,]
betaCovMatrix_p_tot <- modelMatrix[(P+1):nrow(modelMatrix),]

tot_p_alls <- res_p_alls[!duplicated(res_p_alls$gene), colnames(res_p_alls) %in% c('lambda', 'alpha', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')]

res_tot_p <- list(geneAnn = gene_ann, tot = tot_p_alls, beta_snps = betaSnpsMatrix_p_tot, beta_cov = betaCovMatrix_p_tot, 
                    seed = data.frame(type ='tot', value = seed, nfold = nfolds)) 

save(res_tot_p, file = sprintf('%s/resPrior_NotHeritableGenes_%s.RData', outFold, curChrom))






