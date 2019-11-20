#### code written by Lucia Trastulla, e-mail: lucia_trastulla@psych.mpg.de ####
options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(parallel))

parser <- ArgumentParser(description="part1: evaluate res using nested CV + compute alpha/lambda outer fold, without prior")

parser$add_argument("--curChrom", type = "character", help = "chromosome")
parser$add_argument("--genoDat_file", type = "character", help = "genotype data, common path (ending chr specific)")
parser$add_argument("--geneExp_file", type = "character", help = "RNA expression, complete path")
parser$add_argument("--covDat_file", type = "character", help = "covariance file, complete path")
parser$add_argument("--InfoFold", type = "character", help = "path to fold with additional info (e.g. snp-gene dist matrix)")
# parser$add_argument("--priorDat_file", type = "character", help = "prior matrix file, common path (ending chr specific)")
# parser$add_argument("--priorInf", type="integer", nargs = '*' , default = 0, help = "index prior feature to be used, if 0 all column are used")
parser$add_argument("--ncores", type="integer", default = 20, help = "n.of cores for the parallelization")
parser$add_argument("--functR", type="character", help = "Rscript with functions to be used, complete path")
parser$add_argument("--seed_out", type="integer", default = 1234, help = "seed to create outer folders")
parser$add_argument("--seed_in", type="integer", default = 42, help = "seed to create inner folders")
parser$add_argument("--nfolds_out", type="integer", default = 5, help = "n. of outer folds")
parser$add_argument("--nfolds_in", type="integer", default = 5, help = "n. of inner folds")
parser$add_argument("--cis_thres", type = "integer", default = 200000, help = "window (in bp) to compute distance from gene and snps")
parser$add_argument("--Dx", type = "logical", default = FALSE, help = "if true, Dx included as covariate")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
curChrom <- args$curChrom
genoDat_file <- args$genoDat_file
geneExp_file <- args$geneExp_file
covDat_file <- args$covDat_file
ncores <- args$ncores
functR <- args$functR
seed_out <- args$seed_out
seed_in <- args$seed_in
nfolds_out <- args$nfolds_out
nfolds_in <- args$nfolds_in
cis_thres <- args$cis_thres
Dx <- args$Dx
InfoFold <- args$InfoFold
outFold <- args$outFold


# #################################################################
# curChrom <- 'chr21'
# covDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/covariateMatrix_Control50.txt'
# genoDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/Genotyping_data/Genotype_dosage_'
# geneExp_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/RNAseq_data/EXCLUDE_ANCESTRY_SVA/RNAseq_filt.txt'
# ncores <- 5
# seed_out = 1234
# seed_in = 42
# cis_thres = 200000
# nfolds_out = 5
# nfolds_in = 5
# outFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50/'
# InfoFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/'
# functR <- '/ziller/lucia/eQTL_PROJECT_CMC/RSCRIPTS/SCRIPTS_v1/ElNet_withPrior_functions_run.R'
# Dx = F
# #################################################################

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
id_her <- which(gene_ann$type == 'heritable')
expDat <- expDat[id_her, ]
gene_ann <- gene_ann[id_her,]

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
geneSnpDist <- geneSnpDist[, id_her]
covIndex <- c((nrow(geneSnpDist)+1):(nrow(geneSnpDist)+ncol(covDat)))
###################################################################################################################################

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


res <- vector(mode = 'list', length = nfolds_out)
lambdaVec <- vector(mode = 'list', length = nfolds_out)
alphaVec <- vector(mode = 'list', length = nfolds_out)
test_res <- vector(mode = 'list', length = nfolds_out)
betaSnpsMatrix <- vector(mode = 'list', length = nfolds_out)
betaCovMatrix <- vector(mode = 'list', length = nfolds_out)

# not prior info
for(k in 1:nfolds_out){
    
    print(k)
    system.time(res[[k]] <- parLapply(clust, seq(1,N), expPrediction_cv_noPrior_fold, alpha = alpha_h, 
                                      seed = seed_in, nfolds = nfolds_in, fold = setdiff(1:M, Folds[[k]])))
    res[[k]] <- do.call(rbind, res[[k]])
    colnames(res[[k]]) <- c('gene', 'idBeta', 'Beta', 'lambda', 'alpha', 'dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')
    res[[k]] <- as.data.frame(res[[k]])
    print(paste('regression computed fold', k))
    print(head(res[[k]]))
    
    lambdaVec[[k]] <- res[[k]][!duplicated(res[[k]][,1]),4]
    lambdaVec[[k]][is.na(lambdaVec[[k]])] <- 0
    alphaVec[[k]] <- res[[k]][!duplicated(res[[k]][,1]),5]
    alphaVec[[k]][is.na(alphaVec[[k]])] <- 0
    
    # save beta in sparse matrix format
    tmp <- res[[k]][!is.na(res[[k]][,2]),]
    modelMatrix <- sparseMatrix(i=tmp[,2],j=tmp[,1],x=tmp[,3],dims=c(P + ncol(covDat) + 1, N), giveCsparse=T)
    
    betaSnpsMatrix[[k]] <- modelMatrix[1:P,]
    betaCovMatrix[[k]] <- modelMatrix[(P+1):nrow(modelMatrix),]
    
    test_res[[k]] <- deviance_funct(genDat, expDat, covDat, modelMatrix, Folds[[k]])
    colnames(test_res[[k]]) <- c('dev', 'dev_geno','dev_cov', 'dev_geno_cov', 'cor', 'cor_pval')
    test_res[[k]] <- as.data.frame(test_res[[k]])
    
}

predTest_geno_comb <- lapply(1:nfolds_out, function(x) as.matrix(t(betaSnpsMatrix[[k]]) %*% t(genDat[Folds[[x]],])))
predTest_geno_comb <- do.call(cbind, predTest_geno_comb)
predTest_cov_comb <- lapply(1:nfolds_out, function(x) as.matrix(t(betaCovMatrix[[k]]) %*% t(cbind(covDat[Folds[[x]],], rep(1,length(Folds[[x]]))))))

adjExpTest_comb <- lapply(1:nfolds_out, function(x) t(expDat[Folds[[x]],]) - predTest_cov_comb[[x]])
adjExpTest_comb <- do.call(cbind, adjExpTest_comb)

cor_comb_test <- data.frame(cor = rep(0,N), cor_pval = rep(0,N))
for(g in 1:N){
  cor_comb_test$cor[g] <- cor.test(adjExpTest_comb[g,],predTest_geno_comb[g,])$estimate
  cor_comb_test$cor_pval[g] <- cor.test(adjExpTest_comb[g,],predTest_geno_comb[g,])$p.value
}


train_res <- lapply(res, function(x) x[!duplicated(x$gene), colnames(x) %in% c('dev', 'dev_geno', 'dev_cov', 'dev_geno_cov', 'dev_lmgeno', 'cor', 'cor_pval')])

res_tot <- list(geneAnn = gene_ann, train = train_res, test = test_res, cor_comb_test = cor_comb_test,
                beta_snps = betaSnpsMatrix, beta_cov = betaCovMatrix, 
                seed = data.frame(type = c('inner', 'outer'), value = c(seed_in, seed_out), nfold = c(nfolds_in, nfolds_out))) 

lambdaMat <- cbind(gene_ann$ensembl_gene_id, as.data.frame(do.call(cbind, lambdaVec)))
colnames(lambdaMat) <- c('ensembl_gene_id', sapply(1:nfolds_out, function(x)paste0('Fold', x)))

alphaMat <- cbind(gene_ann$ensembl_gene_id, as.data.frame(do.call(cbind, alphaVec)))
colnames(alphaMat) <- c('ensembl_gene_id', sapply(1:nfolds_out, function(x)paste0('Fold', x)))

# store results
write.table(x = lambdaMat, file = sprintf('%s/optim_lambda_%s.txt', outFold,  curChrom), col.names = T, row.names = F, sep = '\t', quote = F)
write.table(x = alphaMat, file = sprintf('%s/optim_alpha_%s.txt', outFold,  curChrom), col.names = T, row.names = F, sep = '\t', quote = F)

save(res_tot, file = sprintf('%s/resNoPrior_NestedCV_HeritableGenes_%s.RData', outFold, curChrom))
  
