#!/usr/bin/env Rscript

#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(nloptr))
suppressPackageStartupMessages(library(bigmemory))


parser <- ArgumentParser(description="part2: find optimal E parameter on the outer fold, use already computed alpha and lambda")

parser$add_argument("--genoDat_file", type = "character", help = "genotype data, common path (ending chr specific)")
parser$add_argument("--geneExp_file", type = "character", help = "RNA expression, complete path")
parser$add_argument("--covDat_file", type = "character", help = "covariance file, complete path")
parser$add_argument("--InfoFold", type = "character", help = "path to fold with additional info (e.g. snp-gene dist matrix)")
parser$add_argument("--part1Res_fold", type = "character", help = "path folder reuslt part1")
parser$add_argument("--priorDat_file", type = "character", help = "prior matrix file, common path (ending chr specific)")
parser$add_argument("--priorInf", type="integer", nargs = '*' , default = 0, help = "index prior feature to be used, if 0 all column are used [default %(default)s]")
parser$add_argument("--ncores", type="integer", default = 10, help = "n.of cores for the parallelization [default %(default)s]")
parser$add_argument("--functR", type="character", help = "Rscript with functions to be used, complete path")
parser$add_argument("--cis_thres", type = "integer", default = 200000, help = "window (in bp) to compute distance from gene and snps [default %(default)s]")
parser$add_argument("--Dx", type = "logical", default = FALSE, help = "if true, Dx included as covariate [default %(default)s]")
parser$add_argument("--maxIter", type="integer", default = 20, help = "maximum number of iterations [default %(default)s]")
parser$add_argument("--dThres", type="double", default = 0.001, help = "difference threshold to stop iteration [default %(default)s]")
parser$add_argument("--E_set", type="double", nargs = '*', help = "range E parameter to search")
parser$add_argument("--convert_par", type="double",  default = 0.25, help = "parameter to rescale all priors, needed to start E par search not too low [default %(default)s]")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
genoDat_file <- args$genoDat_file
geneExp_file <- args$geneExp_file
covDat_file <- args$covDat_file
priorDat_file <- args$priorDat_file
priorInf <- args$priorInf
ncores <- args$ncores
functR <- args$functR
cis_thres <- args$cis_thres
Dx <- args$Dx
maxIter <- args$maxIter
dThres <- args$dThres
E_set <- args$E_set
InfoFold <- args$InfoFold
part1Res_fold <- args$part1Res_fold
convert_par <- args$convert_par
outFold <- args$outFold

# #################################################################
# covDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/covariateMatrix_Control50.txt'
# genoDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/Genotyping_data/Genotype_dosage_'
# geneExp_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/RNAseq_data/EXCLUDE_ANCESTRY_SVA/RNAseq_filt.txt'
# ncores <- 5
# cis_thres = 200000
# outFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/'
# InfoFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/'
# functR <- '/ziller/lucia/eQTL_PROJECT_CMC/RSCRIPTS/SCRIPTS_v1/ElNet_withPrior_functions_run.R'
# Dx = F
# priorDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/priorMatrix_'
# priorInf <- 1:15
# maxIter <- 20
# dThres <- 0.001
# E_set <- 1:5
# part1Res_fold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/200kb/'
# #################################################################

source(functR)
print(paste('Include Dx in covariates:',Dx))

all_Chroms <- paste0('chr', 1:22)
#all_Chroms <- paste0('chr', 1:2)
ngenes <- c()
id_notzero <- list()

# load alpha and lambda file
lambdaMat <- list()
alphaMat <- list()
nvar <- data.frame(chr = all_Chroms, nvar = 0)

for(i in 1:length(all_Chroms)){
  
  chr <- all_Chroms[i]
  id_notzero[[i]] <- data.frame(id = NULL, chr  = NULL)
  
  lambda_file <- sprintf('%soptim_lambda_%s.txt', part1Res_fold, chr)
  print(lambda_file)
  alpha_file <- sprintf('%soptim_alpha_%s.txt', part1Res_fold, chr)
  print(alpha_file)
  
  lambdaVec <- read.table(lambda_file, header = T, stringsAsFactors = F)
  ngenes <- c(ngenes, nrow(lambdaVec))
  id_notzero[[i]] <- data.frame(id = which(lambdaVec$Fold1!=0), chr = chr, stringsAsFactors = F)
  
  lambdaMat[[i]] <- lambdaVec[, !colnames(lambdaVec) %in% 'ensembl_gene_id']
  alphaVec <- read.table(alpha_file, header = T, stringsAsFactors = F)
  alphaMat[[i]] <- alphaVec[, !colnames(alphaVec) %in% 'ensembl_gene_id']
  
  # load seed, nfold_out and nVars
  res_part1 <- get(load(sprintf('%sresNoPrior_NestedCV_HeritableGenes_%s.RData', part1Res_fold, chr)))
  nvar$nvar[i] <- nrow(res_part1$beta_snps[[1]])
  
  
}

nfolds_out <- length(res_part1$beta_snps)
seed_out <-  res_part1$seed$val[res_part1$seed$type == 'outer']

ngenes <- c(0, ngenes)
ngenes <- cumsum(ngenes)

id_notzero <- do.call(c, lapply(1:length(all_Chroms), function(x) id_notzero[[x]]$id + ngenes[x])) 
# id_notzero <- lapply(1:length(all_Chroms), function(x) id_notzero[[x]]$id)
gamma0_it <- rep(0, length(E_set)) 

print(gamma0_it)

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
gene_ann <- expDat[expDat$chrom %in% all_Chroms , colnames(expDat) %in% c('type',	'chrom', 'TSS_start','TSS_end','name','start_position','end_position','ensembl_gene_id','external_gene_name')]
expDat <- expDat[expDat$chrom %in% all_Chroms, !colnames(expDat) %in% c('type',	'chrom', 'TSS_start','TSS_end','name','start_position','end_position','ensembl_gene_id','external_gene_name')]
# order based on sample ann
id <- unname(sapply(sampleAnn$RNASample_ID, function(x) which(x == colnames(expDat))))
expDat <- expDat[, id]
# consider only heritable genes and save the id
id_her <- which(gene_ann$type == 'heritable')
id_her_chr <- sapply(all_Chroms, function(x) which(gene_ann$type[gene_ann$chrom == x] == 'heritable'))
expDat <- expDat[id_her, ]
gene_ann <- gene_ann[id_her,]

## combine info from all chr
genDat <- list()
pDat <- list()
geneSnpDist <- list()

cis_ann <- paste0(cis_thres/10^5, 'e+05')
print(cis_ann)

for(i in 1:length(all_Chroms)){
  
  # genotype
  chr <- all_Chroms[i]
  print(chr)
  
  tmp <- read.table(gzfile(sprintf('%s%s_matrix.txt.gz', genoDat_file, chr)), header = T, sep = '\t', check.names=F)
  # order based on sample ann
  id <- unname(sapply(sampleAnn$genoSample_ID, function(x) which(x == colnames(tmp))))
  tmp <- tmp[, id]
  
  genDat[[i]] <- big.matrix(ncol=nvar$nvar[i] , nrow = ncol(expDat), type = "double", init = 0, dimnames = NULL, shared = F)
  genDat[[i]][,] <- t(as.matrix(tmp))
  
  ## prior matrix
  pDat[[i]] <- read.table(gzfile(sprintf('%s%s.txt.gz', priorDat_file, chr)), header = T, sep = '\t')
  pNames <- colnames(pDat[[i]])[priorInf]
  pDat[[i]] <- as.matrix(pDat[[i]][,priorInf]) # only cell type of interest
  
  ## gene-snp distance matrix
  geneSnpDist[[i]] <- readMM(sprintf('%sENSEMBL_gene_SNP_%s_%s_matrix.mtx', InfoFold, cis_ann, chr))
  geneSnpDist[[i]] <- geneSnpDist[[i]][, id_her_chr[[i]]]
  
}

# normalize epi priors
type_val <- apply(pDat[[1]], 2, unique)
continous_pos <- !sapply(type_val, function(x) all(x %in% c(0,1)) & length(x)==2) 
tmp <- do.call(rbind, pDat)

print(pNames)

for(i in 1:length(all_Chroms)){
  pDat[[i]] <-  pDat[[i]]*convert_par
}
print(head(pDat[[1]]))

tot_var <- sum(sapply(geneSnpDist, nrow))
covIndex_chr <- sapply(nvar$nvar, function(x) c((x+1):(x+ncol(covDat))))
covIndex <- c((tot_var+1):(tot_var+ncol(covDat)))

###################################################################################################################################
expDat <- t(expDat)

##
N <- ncol(expDat) # n. of genes
N_chr <- sapply(all_Chroms, function(x) length(which(gene_ann$chrom == x)))
M <- nrow(expDat) # n. of samples
P_chr <- sapply(genDat,ncol) # n. of SNPs
P <- sum(P_chr)
K <- ncol(pDat[[1]]) # n. of priors
##

# inizialize priors
priorMat <- pDat

########################
######## step 2 ########
########################
# for fixed lambda and alpha optimize weights
# find the best E_h parameter using cross validation

ind_SNPs <- lapply(1:length(all_Chroms), function(y) lapply(1:N_chr[y], function(X) which(geneSnpDist[[y]][,X]!=0)))

# use gradient descent algorithm
set.seed(seed_out)
Folds <- generateCVRuns(nfold = nfolds_out, labels = 1:M, ntimes = 1)[[1]]

expDat_k <- lapply(Folds, function(x) expDat[-x,])
covDat_k <- lapply(Folds, function(x) covDat[-x,])

registerDoParallel(cores=min(ncores, length(E_set)))

# parallelize over E
res_E <- foreach(l=1:length(E_set), .combine='comb', .multicombine=TRUE, 
                 .init=list(list(), list(), list(), list(),list(),list(), list(), list(), list(), list(), list(),list(), list(),
                            list(), list(), list(), list(), list(), list(),list(), list(), list(), list(), list(), list(), list(), list(), list(), list()))%dopar%{
                              
                              E_h <- E_set[l]
                              print(paste0('E parameter:', E_h)) 
                              
                              weightMat <- vector(mode = 'list', length = nfolds_out)
                              fStats_test <- vector(mode = 'list', length = nfolds_out)
                              fStats_train <- vector(mode = 'list', length = nfolds_out)
                              nCountList <-  vector(mode = 'list', length = nfolds_out)
                              devTrain <- vector(mode = 'list', length = nfolds_out)
                              devTrain_geno <- vector(mode = 'list', length = nfolds_out)
                              devTrain_cov <- vector(mode = 'list', length = nfolds_out)
                              devTrain_genocov <- vector(mode = 'list', length = nfolds_out)
                              devTrain_lmgeno <- vector(mode = 'list', length = nfolds_out)
                              corTrain <- vector(mode = 'list', length = nfolds_out)
                              corTrain_pval <- vector(mode = 'list', length = nfolds_out)
                              corTrain_noadj <- vector(mode = 'list', length = nfolds_out)
                              corTrain_noadj_pval <- vector(mode = 'list', length = nfolds_out)
                              SquErr <- vector(mode = 'list', length = nfolds_out)
                              devTest <- vector(mode = 'list', length = nfolds_out)
                              devTest_geno <- vector(mode = 'list', length = nfolds_out)
                              devTest_cov <- vector(mode = 'list', length = nfolds_out)
                              devTest_genocov <- vector(mode = 'list', length = nfolds_out)
                              corTest <- vector(mode = 'list', length = nfolds_out)
                              corTest_pval <- vector(mode = 'list', length = nfolds_out)
                              corTest_noadj <- vector(mode = 'list', length = nfolds_out)
                              corTest_noadj_pval <- vector(mode = 'list', length = nfolds_out)
                              tot_obj <- vector(mode = 'list', length = nfolds_out)
                              tot_obj_test <- vector(mode = 'list', length = nfolds_out)
                              nIter <- rep(0, nfolds_out)
                              predTest_geno_fin <- vector(mode = 'list', length = nfolds_out)
                              predTest_cov_fin <- vector(mode = 'list', length = nfolds_out)
                              betaSnpsMatrix <- vector(mode = 'list', length = nfolds_out)
                              betaCovMatrix <- vector(mode = 'list', length = nfolds_out)
                              
                              for(k in 1:nfolds_out){
                                
                                print(paste('Fold:', k))
                                
                                alphaVec <- lapply(alphaMat, function(x) x[,k])
                                lambdaVec <- lapply(lambdaMat, function(x) x[,k])
                                
                                weightMat[[k]] <- matrix(0,maxIter,K)
                                fStats_test[[k]] <- matrix(0,maxIter,3) 
                                fStats_train[[k]] <- matrix(0,maxIter,3) 
                                nCountList[[k]] <-  matrix(0,maxIter,N) 
                                devTrain[[k]] <- matrix(0,maxIter,N)
                                devTrain_geno[[k]] <- matrix(0,maxIter,N)
                                devTrain_cov[[k]] <- matrix(0,maxIter,N)
                                devTrain_genocov[[k]] <- matrix(0,maxIter,N)
                                devTrain_lmgeno[[k]] <- matrix(0,maxIter,N)
                                corTrain[[k]] <- matrix(0,maxIter,N)
                                corTrain_pval[[k]] <- matrix(0,maxIter,N)
                                corTrain_noadj[[k]] <- matrix(0,maxIter,N)
                                corTrain_noadj_pval[[k]] <- matrix(0,maxIter,N)
                                SquErr[[k]] <- matrix(0,maxIter,N)
                                devTest[[k]] <- matrix(0,maxIter,N)
                                devTest_geno[[k]] <- matrix(0,maxIter,N)
                                devTest_cov[[k]] <- matrix(0,maxIter,N)
                                devTest_genocov[[k]] <- matrix(0,maxIter,N)
                                corTest[[k]] <- matrix(0,maxIter,N)
                                corTest_pval[[k]] <- matrix(0,maxIter,N)
                                corTest_noadj[[k]] <- matrix(0,maxIter,N)
                                corTest_noadj_pval[[k]] <- matrix(0,maxIter,N)
                                # Beta_it <- list()
                                # res_partial <- list()
                              
                                obj <- 0
                                tot_obj[[k]] <- 0
                                tot_obj_test[[k]] <- 0
                                cIter=1
                                converged=F
                                pWeight <- rep(gamma0_it[l], K)
                                print(pWeight)
                                
                                while ((cIter <= maxIter) & (!converged)){
                                  
                                  curPrior  <- lapply(priorMat, function(x) getPrior(x, pWeight))
                                  # clusterExport(clust, "curPrior")
                                  
                                  res_k <- list()
                                  modelMatrix <- list()
                                  for(j in 1:length(all_Chroms)){
                                    
                                    # print(j)
                                    res_k[[j]] <- list()
                                    for(i in 1:N_chr[j]){
                                      
                                      #print(i)
                                      res_k[[j]][[i]] <- expPrediction_chr(expDat_k = expDat_k[[k]][, which(gene_ann$chrom == all_Chroms[j])[i]], 
                                                                           genDat_k = genDat[[j]][-Folds[[k]], ind_SNPs[[j]][[i]]], alpha = alphaVec[[j]][i], lambda = lambdaVec[[j]][i], 
                                                                           covDat_k = covDat_k[[k]], prior = curPrior[[j]][ind_SNPs[[j]][[i]]], d = ind_SNPs[[j]][[i]][1]-1, chr = j, X=i)
                                      
                                    }
                                    
                                    res_k[[j]] <- do.call(rbind, res_k[[j]])
                                    colnames(res_k[[j]]) <- c('gene', 'idBeta', 'Beta', 'lambda', 'alpha',  'SquErr', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval', 'cor_noadj', 'cor_noadj_pval')
                                    res_k[[j]]=res_k[[j]][!is.na(res_k[[j]][,2]),]
                                    
                                    
                                    # store beta results
                                    modelMatrix[[j]] <- sparseMatrix(i=res_k[[j]][,2],j=res_k[[j]][,1],x=res_k[[j]][,3],dims=c(P_chr[j] + ncol(covDat_k[[k]]) + 1, N_chr[j]),giveCsparse=T)
                                    
                                    #####################################
                                    # Beta_it[[cIter]] <- modelMatrix
                                    #####################################
                                    
                                  }
                                  print('regression computed')
                                  
                                  
                                  #print(length(id_notzero))
                                  #print(N)
                                  
                                  devTrain[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'dev']))
                                  devTrain_geno[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'dev_geno']))
                                  devTrain_cov[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'dev_cov']))
                                  devTrain_genocov[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'dev_geno_cov']))
                                  devTrain_lmgeno[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'dev_lmgeno']))
                                  corTrain[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'cor']))
                                  corTrain_pval[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'cor_pval']))
                                  corTrain_noadj[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'cor_noadj']))
                                  corTrain_noadj_pval[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'cor_noadj_pval']))
                                  
                                  SquErr[[k]][cIter,id_notzero] = unlist(sapply(res_k, function(x) x[!duplicated(x[,1]), colnames(x) == 'SquErr']))
                                  
                                  nCountList[[k]][cIter,] = unlist(sapply(1:length(all_Chroms), function(x) sapply(seq(1,N_chr[x]), function(X){return(sum(res_k[[x]][,1]==X,na.rm=T))})))
                                  # delete NA results corresponding to no close SNPs to the gene considered
                                  tmp <- lapply(1:length(all_Chroms), function(j) deviance_funct(genDat[[j]], expDat[,which(gene_ann$chrom == all_Chroms[j])], covDat, modelMatrix[[j]], Folds[[k]]))
                                  devTest[[k]][cIter,] <-  unlist(sapply(tmp, function(x) x[,colnames(x) == 'dev.ratio']))
                                  devTest_geno[[k]][cIter,] <- unlist(sapply(tmp, function(x) x[,colnames(x) == 'dev_geno']))
                                  devTest_cov[[k]][cIter,] <- unlist(sapply(tmp, function(x) x[,colnames(x) == 'dev_cov']))
                                  devTest_genocov[[k]][cIter,] <- unlist(sapply(tmp, function(x) x[,colnames(x) == 'dev_geno_cov']))
                                  corTest[[k]][cIter,] =unlist(sapply(tmp, function(x) x[,colnames(x) == 'cor_est']))
                                  corTest_pval[[k]][cIter,] = unlist(sapply(tmp, function(x) x[,colnames(x) == 'cor_pval']))
                                  corTest_noadj[[k]][cIter,] =unlist(sapply(tmp, function(x) x[,colnames(x) == 'cor_est_noadj']))
                                  corTest_noadj_pval[[k]][cIter,] = unlist(sapply(tmp, function(x) x[,colnames(x) == 'cor_pval_noadj']))
                                  
                                  
                                  #### NOTE: on lirnet the external cordinate descendent algorthim stops only at max iter
                                  # initial weigths are the one previously estimated
                                  res0 <- nloptr(x0 = pWeight, eval_f=getStepFunction_chr, eval_grad_f=getGradient_chr, lb = rep(0, K),
                                                 opts = list("algorithm" = "NLOPT_LD_MMA","xtol_rel"=1.0e-4,"print_level" = 0),
                                                 modelMatrix = modelMatrix, priorMat = priorMat, alphas = alphaVec, lambdas = lambdaVec, E_h = E_h)
                                  
                                  # save the solution as the new weigths
                                  pWeightNew <- res0$solution
                                  pWeight <- pWeightNew
                                  print(pWeight)
                                  
                                  # save weigth
                                  weightMat[[k]][cIter,] <- pWeightNew
                                  
                                  # cv error on the test and train set
                                  names <- c('obj', 'pen', 'pen E_par')
                                  cvErr_test <-  getErrorComponents_chr(lambdaVec, genDat, expDat, covDat, modelMatrix, pWeight, priorMat, alphaVec, E_h, Folds[[k]], all_Chroms, gene_ann)
                                  print(sapply(1:length(names), function(x) paste('test:', names[x], cvErr_test[x])))
                                  
                                  cvErr_train <-  getErrorComponents_chr(lambdaVec, genDat, expDat, covDat, modelMatrix, pWeight, priorMat, alphaVec, E_h, setdiff(1:M, Folds[[k]]), all_Chroms, gene_ann)
                                  print(sapply(1:length(names), function(x) paste('train:', names[x], cvErr_train[x])))
                                  
                                  fStats_test[[k]][cIter,] <- cvErr_test
                                  fStats_train[[k]][cIter,] <- cvErr_train
                                  print(fStats_train[[k]][,1])
                                  print(fStats_test[[k]][,1])
                                  
                                  # update total objective function
                                  objNew <- sum(cvErr_train)
                                  
                                  if ((abs(obj-objNew))<dThres){
                                    converged=T
                                  }
                                  
                                  # update obj
                                  obj <- objNew
                                  tot_obj[[k]] <- c(tot_obj[[k]], obj)
                                  # print(tot_obj)
                                  tot_obj_test[[k]] <- c(tot_obj_test[[k]], sum(cvErr_test))
                                  # print(tot_obj_test)
                                  # update iter
                                  cIter=cIter+1
                                  
                                  # save prediction on test for cor combined, (only for the final values)
                                  if(!((cIter <= maxIter) & (!converged))){
                                    
                                    betaSnpsMatrix[[k]] <- mapply(function(x, y) x[1:y, ], x = modelMatrix, y = P_chr)
                                    betaCovMatrix[[k]] <-  mapply(function(x, y) x[(y+1):nrow(x), ], x = modelMatrix, y = P_chr)
                                    predTest_geno_fin[[k]] <- lapply(1:length(all_Chroms), function(x) as.matrix(t(betaSnpsMatrix[[k]][[x]]) %*% t(genDat[[x]][Folds[[k]],])))
                                    predTest_geno_fin[[k]] <- do.call(rbind, predTest_geno_fin[[k]])
                                    predTest_cov_fin[[k]] <- lapply(1:length(all_Chroms), function(x) as.matrix(t(betaCovMatrix[[k]][[x]]) %*% t(cbind(covDat[Folds[[k]],], rep(1,length(Folds[[k]]))))))
                                    predTest_cov_fin[[k]] <- do.call(rbind, predTest_cov_fin[[k]])
                                    
                                  }
                                  
                                  
                                }
                                
                                nIter[k] <- cIter - 1
                                
                              }
                              
                              list(weightMat, fStats_train, fStats_test, nCountList, SquErr, tot_obj, tot_obj_test, nIter, devTrain, devTrain_geno, 
                                   devTrain_cov, devTrain_genocov, devTrain_lmgeno, corTrain, corTrain_pval, corTrain_noadj, corTrain_noadj_pval, devTest, devTest_geno, devTest_cov,
                                   devTest_genocov, corTest, corTest_pval, corTest_noadj, corTest_noadj_pval, predTest_geno_fin, predTest_cov_fin, betaSnpsMatrix, betaCovMatrix)
                              
                            }  

# save total result
save(res_E, file = sprintf('%sresE_allchr.RData', outFold))

## trasnform output in the correct format
nIter <- res_E[[8]]

obj_steps_train <- lapply(res_E[[6]], function(x) lapply(x, function(y) y[-1]))
obj_steps_test <- lapply(res_E[[7]], function(x) lapply(x, function(y) y[-1]))

evalf_steps_train <- lapply(res_E[[2]], function(x) lapply(x, function(y) y[,1]))
evalf_steps_train <- lapply(1:length(E_set), function(l) lapply(1:nfolds_out, function(k) evalf_steps_train[[l]][[k]][1:nIter[[l]][k]]))

evalf_steps_test <- lapply(res_E[[3]], function(x) lapply(x, function(y) y[,1]))
evalf_steps_test <- lapply(1:length(E_set), function(l) lapply(1:nfolds_out, function(k) evalf_steps_test[[l]][[k]][1:nIter[[l]][k]]))

#####
cv_train_err <- matrix(0, nrow = nfolds_out, ncol = length(E_set))
cv_test_err <- matrix(0, nrow = nfolds_out, ncol = length(E_set))

dev_train <- vector(mode = 'list', length = length(E_set))
dev_geno_train <- vector(mode = 'list', length = length(E_set))
dev_cov_train <- vector(mode = 'list', length = length(E_set))
dev_geno_cov_train <- vector(mode = 'list', length = length(E_set))
dev_lmgeno_train <- vector(mode = 'list', length = length(E_set))
cor_train <-  vector(mode = 'list', length = length(E_set))
cor_pval_train <-  vector(mode = 'list', length = length(E_set))
cor_noadj_train <-  vector(mode = 'list', length = length(E_set))
cor_noadj_pval_train <-  vector(mode = 'list', length = length(E_set))

dev_test <- vector(mode = 'list', length = length(E_set))
dev_geno_test <- vector(mode = 'list', length = length(E_set))
dev_cov_test <- vector(mode = 'list', length = length(E_set))
dev_geno_cov_test <- vector(mode = 'list', length = length(E_set))
cor_test <- vector(mode = 'list', length = length(E_set))
cor_pval_test <- vector(mode = 'list', length = length(E_set))
cor_noadj_test <- vector(mode = 'list', length = length(E_set))
cor_noadj_pval_test <- vector(mode = 'list', length = length(E_set))

cor_test_comb <-  vector(mode = 'list', length = length(E_set))
cor_test_noadj_comb <-  vector(mode = 'list', length = length(E_set))

betaSnpsMatrix <- res_E[[28]]
betaCovMatrix <- res_E[[29]]


weights_res <-    vector(mode = 'list', length = length(E_set))

for(l in 1:length(E_set)){
  
  weights_res[[l]] <- sapply(1:nfolds_out, function(k) res_E[[1]][[l]][[k]][nIter[[l]][k],])
  
  cv_train_err[,l] <- sapply(1:nfolds_out, function(k) evalf_steps_train[[l]][[k]][nIter[[l]][k]])
  
  dev_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[9]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  dev_geno_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[10]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  dev_cov_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[11]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  dev_geno_cov_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[12]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  dev_lmgeno_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[13]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  
  cor_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[14]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  cor_pval_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[15]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  cor_noadj_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[16]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  cor_noadj_pval_train[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[17]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  
  
  cv_test_err[,l] <- sapply(1:nfolds_out, function(k) evalf_steps_test[[l]][[k]][nIter[[l]][k]])
  
  dev_test[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[18]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  dev_geno_test[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[19]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  dev_cov_test[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[20]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  dev_geno_cov_test[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[21]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  cor_test[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[22]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  cor_pval_test[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[23]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  cor_noadj_test[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[24]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  cor_noadj_pval_test[[l]] <- do.call(cbind, lapply(1:nfolds_out, function(k) res_E[[25]][[l]][[k]][nIter[[l]][k],])) # dim Nxnfolds
  
  pred_test <- do.call(cbind, res_E[[26]][[l]])
  adj_exp <- lapply(1:nfolds_out, function(x) t(expDat[Folds[[x]],]) - res_E[[27]][[l]][[x]])
  adj_exp <- do.call(cbind, adj_exp)
  exp <- lapply(1:nfolds_out, function(x) t(expDat[Folds[[x]],]))
  exp <- do.call(cbind, exp)
  
  cor_test_comb[[l]] <- data.frame(cor = rep(0,N), cor_pval = rep(0,N))
  cor_test_noadj_comb[[l]] <- data.frame(cor = rep(0,N), cor_pval = rep(0,N))
  for(g in 1:N){
    cor_test_comb[[l]]$cor[g] <- cor.test(adj_exp[g,],pred_test[g,])$estimate
    cor_test_comb[[l]]$cor_pval[g] <- cor.test(adj_exp[g,],pred_test[g,])$p.value
    
    cor_test_noadj_comb[[l]]$cor[g] <- cor.test(exp[g,],pred_test[g,])$estimate
    cor_test_noadj_comb[[l]]$cor_pval[g] <- cor.test(exp[g,],pred_test[g,])$p.value
  }
  
}

colnames(cv_train_err) <- sapply(E_set, function(x) sprintf('E_h:%.2f', x))
colnames(cv_test_err) <- sapply(E_set, function(x) sprintf('E_h:%.2f', x))

write.table(x = cv_train_err, file = sprintf('%scv_train_Epar_allchr.txt', outFold), sep = '\t', row.names = F, col.names = T, quote = F)
write.table(x = cv_test_err, file = sprintf('%scv_test_Epar_allchr.txt', outFold), sep = '\t', row.names = F, col.names = T, quote = F)

# save tot_obj (for all the possible parameters, check the obj is decreasing)
names(obj_steps_train) <- sapply(E_set, function(x) sprintf('E_h:%.2f', x))
names(obj_steps_test) <- sapply(E_set, function(x) sprintf('E_h:%.2f', x))
names(evalf_steps_train) <- sapply(E_set, function(x) sprintf('E_h:%.2f', x))
names(evalf_steps_test) <- sapply(E_set, function(x) sprintf('E_h:%.2f', x))

save(obj_steps_test, file = sprintf('%sobj_Epar_cvtest_allchr.RData', outFold))
save(obj_steps_train, file = sprintf('%sobj_Epar_cvtrain_allchr.RData', outFold))  
save(evalf_steps_test, file = sprintf('%sevalf_Epar_cvtest_allchr.RData', outFold))
save(evalf_steps_train, file = sprintf('%sevalf_Epar_cvtrain_allchr.RData', outFold))  

# optimal E parameter
if(length(E_set)>1){
  id_opt <- which.min(colMeans(cv_test_err))  
}else{
  id_opt <- 1
}

# if the minimum is reached in the interval, use that otherwise
# use as optimal value the one that lead to convergence: first value st diff lower than 0.5
if(length(E_set)>1){
  if(which.min(colMeans(cv_test_err))<ncol(cv_test_err)){
    id_opt <- which.min(colMeans(cv_test_err))
    E_hat <- as.numeric(strsplit(names(id_opt), 'E_h.')[[1]][2])
    E_par <- -1
    print(sprintf('optimal E parameter: %.2f at position %i', E_hat, id_opt))
  }else{
    id_opt <- which(abs(diff(colMeans(cv_test_err))) < 0.5)[1] + 1
    E_hat <- as.numeric(strsplit(names(id_opt), 'E_h.')[[1]][2])
    E_par <- E_hat
    print(sprintf('minimum not reached in the interval, E parameter based on convergence: %.2f at position %i', E_par, id_opt))
  }
}else{
  id_opt <- 1
  E_hat <- E_set
  E_par <- E_hat
}


weights_opt <-  matrix(weights_res[[id_opt]], ncol=nfolds_out, nrow = length(pNames))
rownames(weights_opt) <- pNames
colnames(weights_opt) <-  sapply(1:nfolds_out, function(x) sprintf('Fold_%i', x))

train_res_opt <- lapply(1:nfolds_out, function(x) data.frame(dev = dev_train[[id_opt]][,x], dev_geno = dev_geno_train[[id_opt]][,x], 
                                                             dev_cov = dev_cov_train[[id_opt]][,x], dev_geno_cov = dev_geno_cov_train[[id_opt]][,x], 
                                                             dev_lmgeno = dev_lmgeno_train[[id_opt]][,x], 
                                                             cor = cor_train[[id_opt]][,x],cor_pval = cor_pval_train[[id_opt]][,x], 
                                                             cor_noadj = cor_noadj_train[[id_opt]][,x], cor_noadj_pval = cor_noadj_pval_train[[id_opt]][,x]))

test_res_opt <- lapply(1:nfolds_out, function(x) data.frame(dev = dev_test[[id_opt]][,x], dev_geno = dev_geno_test[[id_opt]][,x], 
                                                            dev_cov = dev_cov_test[[id_opt]][,x], dev_geno_cov = dev_geno_cov_test[[id_opt]][,x], 
                                                            cor = cor_test[[id_opt]][,x], cor_pval = cor_pval_test[[id_opt]][,x], 
                                                            cor_noadj = cor_noadj_test[[id_opt]][,x], cor_noadj_pval = cor_noadj_pval_test[[id_opt]][,x]))

res_tot_opt <- list(geneAnn = gene_ann, train_opt = train_res_opt, test_opt = test_res_opt, 
                    cor_comb_test_opt = cor_test_comb[[id_opt]], cor_comb_test_noadj_opt = cor_test_noadj_comb[[id_opt]],
                    beta_snps_opt = betaSnpsMatrix[[id_opt]], beta_cov_opt = betaCovMatrix[[id_opt]], weights_opt = weights_opt, E_opt = E_hat) 


if(E_par<0){
  save(x = res_tot_opt, file = sprintf('%sresPrior_EOpt_NestedCV_HeritableGenes_allchr.RData', outFold))  
}else{
  save(x = res_tot_opt, file = sprintf('%sresPrior_EFixed%.2f_NestedCV_HeritableGenes_allchr.RData', outFold, E_par))
}




