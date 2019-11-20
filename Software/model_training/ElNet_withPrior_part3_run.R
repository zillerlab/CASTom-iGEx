#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(bigmemory))
suppressPackageStartupMessages(library(nloptr))



parser <- ArgumentParser(description="part3: With optimal E parameter, create final models with (and without) prior")

parser$add_argument("--genoDat_file", type = "character", help = "genotype data, common path (ending chr specific)")
parser$add_argument("--geneExp_file", type = "character", help = "RNA expression, complete path")
parser$add_argument("--covDat_file", type = "character", help = "covariance file, complete path")
parser$add_argument("--InfoFold", type = "character", help = "path to fold with additional info (e.g. snp-gene dist matrix)")
parser$add_argument("--part2Res_fold", type = "character", help = "path folder reuslt part2")
parser$add_argument("--priorDat_file", type = "character", help = "prior matrix file, common path (ending chr specific)")
parser$add_argument("--priorInf", type="integer", nargs = '*' , default = 0, help = "index prior feature to be used, if 0 all column are used")
parser$add_argument("--ncores", type="integer", default = 10, help = "n.of cores for the parallelization")
parser$add_argument("--functR", type="character", help = "Rscript with functions to be used, complete path")
parser$add_argument("--cis_thres", type = "integer", default = 200000, help = "window (in bp) to compute distance from gene and snps")
parser$add_argument("--Dx", type = "logical", default = FALSE, help = "if true, Dx included as covariate")
parser$add_argument("--maxIter", type="integer", default = 20, help = "maximum number of iterations")
parser$add_argument("--dThres", type="double", default = 0.001, help = "difference threshold to stop iteration")
parser$add_argument("--seed", type="integer", default = 4321, help = "seed to create folders")
parser$add_argument("--nfolds", type="integer", default = 5, help = "n. of folds")
parser$add_argument("--convert_par", type="double",  default = 0.25, help = "")
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
InfoFold <- args$InfoFold
part2Res_fold <- args$part2Res_fold
seed <- args$seed
nfolds <- args$nfolds
convert_par <- args$convert_par
outFold <- args$outFold

# #################################################################
# covDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/covariateMatrix_Control50.txt'
# genoDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/Genotyping_data/Genotype_dosage_'
# geneExp_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/RNAseq_data/EXCLUDE_ANCESTRY_SVA/RNAseq_filt.txt'
# ncores <- 10
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
# part2Res_fold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/200kb/'
# nfolds <- 5
# seed <- 4321
# #################################################################


##################################
### load best E par parameters ###
##################################

source(functR)
all_Chroms <- paste0('chr', 1:22)

# find optimal E parameter
cv_test_file <-sprintf('%s/cv_test_Epar_allchr.txt', part2Res_fold)
cv_test <- read.table(cv_test_file, header = T, stringsAsFactors = F)

# if the minimum is reached in the interval, use that otherwise
# use as optimal value the one that lead to convergence: first value st diff lower than 0.5
if(which.min(colMeans(cv_test))<ncol(cv_test)){
  id_opt <- which.min(colMeans(cv_test))
  E_hat <- as.numeric(strsplit(names(id_opt), 'E_h.')[[1]][2])
  E_par = -1
}else{
  id_opt <- which(abs(diff(colMeans(cv_test))) < 0.5)[1] + 1
  E_hat <- as.numeric(strsplit(names(id_opt), 'E_h.')[[1]][2])
  E_par <- E_hat
}


print(E_hat)
print(id_opt)
print(E_par)

# load nvar
resOpt <- get(load(sprintf('%s/resPrior_EOpt_NestedCV_HeritableGenes_allchr.RData', part2Res_fold)))
nvar <- data.frame(chr = all_Chroms, nvar = sapply(resOpt$beta_snps_opt[[1]],  nrow))


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
  if(length(priorInf) == 1 & priorInf[1] == 0){priorInf <- 1:ncol(pDat[[i]])}
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
# pNames <- colnames(tmp)

print(pNames)

#if(length(which(!continous_pos))>1){
#  N_k <- apply(tmp[,!continous_pos], 2, function(x) length(which(x != 0)))
#  pDat <- lapply(pDat, function(x) cbind(x[,continous_pos], sweep(x[,!continous_pos], 2, N_k/min(N_k), '/')))
#  # pDat[, !gwas_pos] <- sweep(pDat[,!gwas_pos], 2, N_k/min(N_k), '/')
#}

# N_k <- apply(tmp, 2, function(x) mean(x))
# pDat <- lapply(pDat, function(x)  sweep(x, 2, N_k/min(N_k), '/'))
for(i in 1:length(all_Chroms)){
  pDat[[i]] <-  pDat[[i]]*convert_par
}
print(head(pDat[[1]]))


tot_var <- sum(sapply(geneSnpDist, nrow))
covIndex_chr <- sapply(nvar$nvar, function(x) c((x+1):(x+ncol(covDat))))
covIndex <- c((tot_var+1):(tot_var+ncol(covDat)))

###################################################################################################################################
# input values
alpha_h <- seq(0.1,0.9,0.1) # exclude lasso because it choses only 1 varaible different from zero between the correlated ones 

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

#############################################################
### compute lambda/alpha without prior on the entire set ####
#############################################################

# use iterative coordinate descent
clust <- makeCluster(ncores)
clusterExport(clust, "geneSnpDist")
clusterExport(clust, "covIndex_chr")
clusterExport(clust, "covDat")
clusterExport(clust, "expDat")
clusterExport(clust, "genDat")
clusterExport(clust, "gene_ann")
clusterExport(clust, "all_Chroms")
clusterEvalQ(clust,library(glmnet))
clusterExport(clust,'generateCVRuns')

res_noprior <- vector(mode = 'list', length = length(all_Chroms))
lambdaVec <- vector(mode = 'list', length = length(all_Chroms))
alphaVec <- vector(mode = 'list', length = length(all_Chroms))
betaSnpsMatrix_noprior <- vector(mode = 'list', length = length(all_Chroms))
betaCovMatrix_noprior <- vector(mode = 'list', length = length(all_Chroms))
ngenes <-  c()
id_notzero <-  vector(mode = 'list', length = length(all_Chroms))

# not prior info
for(chr in 1:length(all_Chroms)){
  
  print(chr)
  system.time(res_noprior[[chr]] <- parLapply(clust, seq(1,N_chr[chr]), expPrediction_cv_noPrior_chr, alpha = alpha_h, chr = chr, genDat = genDat[[chr]][,],
                                    seed = seed, nfolds = nfolds))
  
  res_noprior[[chr]]  <- do.call(rbind, res_noprior[[chr]] )
  colnames(res_noprior[[chr]] ) <- c('gene', 'idBeta', 'Beta', 'lambda', 'alpha', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')
  res_noprior[[chr]]  <- as.data.frame(res_noprior[[chr]] )
  print(paste0('regression computed chr', chr))
  
  lambdaVec[[chr]] <- res_noprior[[chr]][!duplicated(res_noprior[[chr]][,1]),4]
  lambdaVec[[chr]][is.na(lambdaVec[[chr]])] <- 0
  alphaVec[[chr]] <- res_noprior[[chr]][!duplicated(res_noprior[[chr]][,1]),5]
  alphaVec[[chr]][is.na(alphaVec[[chr]])] <- 0
  
  ngenes <- c(ngenes, length(lambdaVec[[chr]]))
  id_notzero[[chr]] <- data.frame(id = which(lambdaVec[[chr]]!=0), chr = paste0('chr', chr), stringsAsFactors = F)
  
  
  # save beta in sparse matrix format
  tmp <- res_noprior[[chr]][!is.na(res_noprior[[chr]][,2]),]
  modelMatrix <- sparseMatrix(i=tmp[,2],j=tmp[,1],x=tmp[,3],dims=c(P_chr[chr] + ncol(covDat) + 1, N_chr[chr]), giveCsparse=T)
  
  betaSnpsMatrix_noprior[[chr]] <- modelMatrix[1:P_chr[chr],]
  betaCovMatrix_noprior[[chr]] <- modelMatrix[(P_chr[chr]+1):nrow(modelMatrix),]
  

}

ngenes <- c(0, ngenes)
ngenes <- cumsum(ngenes)

id_notzero <- do.call(c, lapply(1:length(all_Chroms), function(x) id_notzero[[x]]$id + ngenes[x]))

allsamples_res <- lapply(res_noprior, function(x) x[!duplicated(x$gene), colnames(x) %in% c('lambda', 'alpha', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')])
allsamples_res <- do.call(rbind, allsamples_res) # collapse all chr together

res_tot_nop <- list(geneAnn = gene_ann, tot = allsamples_res, beta_snps = betaSnpsMatrix_noprior, beta_cov = betaCovMatrix_noprior, 
                    seed = data.frame(type ='tot', value = seed, nfold = nfolds)) 

save(res_tot_nop, file = sprintf('%s/resNoPrior_HeritableGenes_allchr.RData', outFold))

##################
### final step ###
##################
# use the optimal E_h parameter found in step 2,
# use the optimal lambda and alpha found before,
# iterate to find the optimal beta and pWeigth on the entire dataset


weightMat <- matrix(0,maxIter,K)
fStats <- matrix(0,maxIter,3) 
nCountList <-  matrix(0,maxIter,N) 
devFin <- matrix(0,maxIter,N) 
devFin_geno <- matrix(0,maxIter,N) 
devFin_cov <- matrix(0,maxIter,N) 
devFin_genocov <- matrix(0,maxIter,N) 
devFin_lmgeno <- matrix(0,maxIter,N)
corFin <- matrix(0,maxIter,N)
corFin_pval <- matrix(0,maxIter,N)
SquErr <- matrix(0,maxIter,N)
BetaIt <- vector(mode = 'list', length = maxIter)

obj <- 0
tot_obj <- NULL
cIter=1
converged=F
gamma0_it <- 0
pWeight <- rep(gamma0_it, K)

clusterExport(clust, "lambdaVec")
clusterExport(clust, "alphaVec")

while ((cIter <= maxIter) & (!converged)){
  
  curPrior  <- lapply(priorMat, function(x) getPrior(x, pWeight))
  
  res <- list()
  modelMatrix <- list()
  for(j in 1:length(all_Chroms)){
    
    print(j)
    # system.time(res[[j]] <- lapply(seq(1,N_chr[j]), function(x) expPrediction_fin_chr(X = x, prior = curPrior[[j]], id_chr = j))) # not parellalized
    system.time(res[[j]] <- parLapply(clust, seq(1,N_chr[j]), expPrediction_fin_chr, prior = curPrior[[j]], id_chr = j, genDat = genDat[[j]][,]))
    
    res[[j]] <- do.call(rbind, res[[j]])
    colnames( res[[j]] ) <- c('gene', 'idBeta', 'Beta', 'lambda', 'alpha',  'SquErr', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')
    res[[j]]=res[[j]][!is.na(res[[j]][,2]),]
    
    # store beta results
    modelMatrix[[j]] <- sparseMatrix(i=res[[j]][,2],j=res[[j]][,1],x=res[[j]][,3],dims=c(P_chr[j] + ncol(covDat) + 1, N_chr[j]),giveCsparse=T)
    
    
  }
  print('regression computed')
  
  SquErr[cIter, id_notzero] = unlist(sapply(res, function(x) x[!duplicated(x[,1]), colnames(x) == 'SquErr']))
  
  devFin[cIter,id_notzero] = unlist(sapply(res, function(x) x[!duplicated(x[,1]),  colnames(x) == 'dev']))
  devFin_geno[cIter,id_notzero] = unlist(sapply(res, function(x) x[!duplicated(x[,1]),  colnames(x) == 'dev_geno']))
  devFin_cov[cIter, id_notzero] = unlist(sapply(res, function(x) x[!duplicated(x[,1]),  colnames(x) == 'dev_cov']))
  devFin_genocov[cIter,id_notzero] = unlist(sapply(res, function(x) x[!duplicated(x[,1]),  colnames(x) == 'dev_geno_cov']))
  devFin_lmgeno[cIter,id_notzero] = unlist(sapply(res, function(x) x[!duplicated(x[,1]),  colnames(x) == 'dev_lmgeno']))
  
  corFin[cIter,id_notzero] <- unlist(sapply(res, function(x) x[!duplicated(x[,1]),  colnames(x) == 'cor']))
  corFin_pval[cIter,id_notzero] <- unlist(sapply(res, function(x) x[!duplicated(x[,1]),  colnames(x) == 'cor_pval']))
  
  nCountList[cIter,] = unlist(sapply(1:length(all_Chroms), function(x) sapply(seq(1,N_chr[x]), function(X){return(sum(res[[x]][,1]==X,na.rm=T))})))
  # delete NA results corresponding to no close SNPs to the gene considered
  
  BetaIt[[cIter]] <- modelMatrix
  
  #### NOTE: on lirnet the external cordinate descendent algorthim stops only at max iter
  # initial weigths are the one previously estimated
  res0 <- nloptr(x0 = pWeight, eval_f=getStepFunction_chr, eval_grad_f=getGradient_chr, lb = rep(0, K),
                 opts = list("algorithm" = "NLOPT_LD_MMA","xtol_rel"=1.0e-4,"print_level" = 3),
                 modelMatrix = modelMatrix, priorMat = priorMat, alphas = alphaVec, lambdas = lambdaVec, E_h = E_hat)
  
  # save the solution as the new weigths
  pWeightNew <- res0$solution
  pWeight <- pWeightNew
  
  # save weigth
  weightMat[cIter,] <- pWeightNew
  
  Err <-  getErrorComponents_chr(lambdaVec, genDat, expDat, covDat, modelMatrix, pWeight, priorMat, alphaVec, E_hat, 1:M, all_Chroms, gene_ann)
  fStats[cIter,] <- Err
  
  # update total objective function
  objNew <- sum(Err)
  
  if ((abs(obj-objNew))<dThres){
    converged=T
  }
  
  # update obj
  obj <- objNew
  tot_obj <- c(tot_obj, obj)
  
  # update iter
  cIter=cIter+1
  
}

nIter <- cIter - 1

################################
###### save final results ######
################################

# iteration steps
weightMat <- weightMat[1:nIter,]
fStats <- fStats[1:nIter,]
nCountList <- nCountList[1:nIter, ]

devFin <- devFin[1:nIter, ]
devFin_geno <- devFin_geno[1:nIter, ]
devFin_cov <- devFin_cov[1:nIter, ]
devFin_genocov <- devFin_genocov[1:nIter, ]
devFin_lmgeno <- devFin_lmgeno[1:nIter, ]
corFin <- corFin[1:nIter, ]
corFin_pval <- corFin_pval[1:nIter,]

SquErr <- SquErr[1:nIter, ]
BetaIt <- lapply(1:nIter, function(x) BetaIt[[x]])


results_iteration <- list(Eopt = E_hat, pWeight = weightMat, errComp = fStats, nCount = nCountList, 
                          dev = devFin, dev_geno = devFin_geno, dev_cov = devFin_cov, dev_genocov = devFin_genocov, cor = corFin, cor_pval = corFin_pval, 
                          MSE = SquErr, beta = BetaIt, obj = tot_obj)

if(E_par<0){
  save(x = results_iteration, file = sprintf('%s/resPrior_EOpt_Iteration_HeritableGenes_allchr.RData', outFold))  
}else{
  save(x = results_iteration, file = sprintf('%s/resPrior_EFixed%.2f_Iteration_HeritableGenes_allchr.RData', outFold, E_par))
}

# final result
betaSnpsMatrix <-  lapply(1:length(all_Chroms), function(x) modelMatrix[[x]][1:P_chr[x],])
betaCovMatrix <-  lapply(1:length(all_Chroms), function(x) modelMatrix[[x]][(P_chr[x]+1):nrow(modelMatrix[[x]]),])
weights_fin <- pWeight
prior_coeff <- lapply(priorMat, function(x) getPrior(x, weights_fin))
tot_p <- data.frame(dev = devFin[nIter,], dev_geno =  devFin_geno[nIter,], dev_cov = devFin_cov[nIter,], dev_genocov = devFin_genocov[nIter, ], dev_lmgeno = devFin_lmgeno[nIter, ], 
                    cor = corFin[nIter, ], cor_pval = corFin_pval[nIter, ])


res_tot_p <- list(geneAnn = gene_ann, tot = tot_p, beta_snps = betaSnpsMatrix, beta_cov = betaCovMatrix, 
                  Eopt = E_hat, weights = weights_fin, prior_coeff = prior_coeff)
  

if(E_par<0){
  save(x = res_tot_p, file = sprintf('%s/resPrior_EOpt_HeritableGenes_allchr.RData', outFold))  
}else{
  save(x = res_tot_p, file = sprintf('%s/resPrior_EFixed%.2f_HeritableGenes_allchr.RData', outFold, E_par))
}



