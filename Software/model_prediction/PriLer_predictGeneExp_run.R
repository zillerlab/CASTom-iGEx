#!/usr/bin/env Rscript

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(bigmemory))

parser <- ArgumentParser(description="Predict gene expression on a new dataset (all genes)")

parser$add_argument("--outTrain_fold", type = "character", help = "path output model")
parser$add_argument("--genoDat_file", type = "character", help = "directory with genotype data to use for the prediction, same format as input for train")
parser$add_argument("--covDat_file", type = "character", help = "directory with sample info for new genotype data, keep only this samples")
parser$add_argument("--InfoFold", type = "character", help = "path to fold with additional info (e.g. snp-gene dist matrix)")
parser$add_argument("--cis_thres", type = "integer", default = 200000, help = "window (in bp) to compute distance from gene and snps [default %(default)s]")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
outTrain_fold <- args$outTrain_fold
genoDat_file <- args$genoDat_file
covDat_file <- args$covDat_file
cis_thres <- args$cis_thres
InfoFold <- args$InfoFold
outFold <- args$outFold

#####################################################################
# outFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/predict_All_oldScale/train_Control50/200kb/'
# outTrain_fold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/train_Control50_oldScale/200kb/'
# genoDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/Genotyping_data/Genotype_dosage_'
# covDat_file <- '/ziller/lucia/eQTL_PROJECT_CMC/INPUT_DATA_SCRIPTS_v1/Covariates/covariateMatrix.txt'
# functR <- '/ziller/lucia/eQTL_PROJECT_CMC/RSCRIPTS/SCRIPTS_v1/ElNet_withPrior_functions_run.R'
# cis_thres <- 200000
# InfoFold <- '/ziller/lucia/eQTL_PROJECT_CMC/OUTPUT_CMC_SCRIPTS_v1/'
#####################################################################

cis_ann <- paste0(cis_thres/10^5, 'e+05')

covDat <- read.table(covDat_file, header = T, sep = '\t', stringsAsFactors = F, check.names=F)
sampleAnn <- covDat[, colnames(covDat) %in% c('Individual_ID', 'genoSample_ID')]
M <- nrow(sampleAnn)

#### load model ####
all_Chroms <- paste0('chr',1:22)
geneInfo <- read.table(sprintf('%sresPrior_regEval_allchr.txt', outTrain_fold), header = T, stringsAsFactors = F)
betaCoeff <- get(load(sprintf('%sresPrior_regCoeffSnps_allchr.RData', outTrain_fold)))


# predict for each chromosome
pExp <- vector(mode = 'list', length = length(all_Chroms))

for(i in 1:length(all_Chroms)){
  
  chr <- all_Chroms[i]
  print(chr)
  
  genDat <- read.table(gzfile(sprintf('%s%s_matrix.txt.gz', genoDat_file, chr)), header  = T, check.names=F)
  # order based the sample_id file (NOTE: covDat already ordered based on that)
  id <- unname(sapply(sampleAnn$genoSample_ID, function(x) which(x == colnames(genDat))))
  genDat <- genDat[, id]
  
  geneSnpDist <- readMM(sprintf('%sENSEMBL_gene_SNP_%s_%s_matrix.mtx', InfoFold, cis_ann, chr))
  ind_SNPs <- lapply(1:ncol(geneSnpDist), function(X) which(geneSnpDist[,X]!=0))
  beta_chr <- betaCoeff[[i]]
  
  N_chr <- ncol(beta_chr)
  pExp[[i]] <- matrix(nrow = N_chr, ncol = M)
  
  for(n in 1:N_chr){
    
    print(n)
    tmp <- genDat[ind_SNPs[[n]], ]
    tmp_beta <- beta_chr[ind_SNPs[[n]], n]
    pExp[[i]][n,] <- tmp_beta %*% as.matrix(tmp)
    
  }
  
  rm(genDat)
  
}

pExp <- do.call(rbind, pExp)
colnames(pExp) <- sampleAnn$Individual_ID
pExp_fin <- cbind(geneInfo, pExp)

# save
write.table(file = sprintf('%spredictedExpression.txt', outFold), x = pExp_fin, col.names = T, row.names = F, sep  = '\t', quote  = F)
system(paste("gzip",sprintf('%spredictedExpression.txt', outFold)))




