#!/usr/bin/env Rscript

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(bigmemory))
suppressPackageStartupMessages(library(data.table))

parser <- ArgumentParser(description="Predict gene expression on a new dataset (all genes)")

parser$add_argument("--outTrain_fold", type = "character", help = "path output model")
parser$add_argument("--genoDat_file", type = "character", help = "directory with genotype data to use for the prediction, same format as input for train")
parser$add_argument("--covDat_file", type = "character", help = "directory with sample info for new genotype data, keep only this samples")
parser$add_argument("--InfoFold", type = "character", help = "path to fold with additional info (e.g. snp-gene dist matrix)")
parser$add_argument("--cis_thres", type = "integer", default = 200000, help = "window (in bp) to compute distance from gene and snps [default %(default)s]")
parser$add_argument("--no_zip", type="logical",default = F, help = "if true final results not zipped [default %(default)s]")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")
parser$add_argument("--genoInfo_model_file", type="character", help = "directory with genotype info data to use for the prediction")
parser$add_argument("--genoInfo_file", type="character", help = "directory with genotype info data from the model, MUST match the trained model")

args <- parser$parse_args()
outTrain_fold <- args$outTrain_fold
genoDat_file <- args$genoDat_file
covDat_file <- args$covDat_file
cis_thres <- args$cis_thres
InfoFold <- args$InfoFold
no_zip <- args$no_zip
outFold <- args$outFold
genoInfo_model_file <- args$genoInfo_model_file
genoInfo_file <- args$genoInfo_file

#####################################################################
# outFold <- './'
# outTrain_fold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/Liver/200kb/CAD_GWAS_bin5e-2/'
# genoDat_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Genotyping_data/German1/Genotype_dosage_'
# covDat_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates//German1/covariateMatrix.txt'
# cis_thres <- 200000
# InfoFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/train_GTEx/Liver/'
# genoInfo_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Genotyping_data/German1/G1.Genotype_VariantsInfo_matchedCADall-UKBB-GTEx_'
# genoInfo_model_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Genotyping_data/UKBB/UKBB.Genotype_VariantsInfo_matchedCADall-UKBB-GTEx_'
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
  
  genDat <- fread(sprintf('%s%s_matrix.txt.gz', genoDat_file, chr), header  = T, check.names=F, data.table = F)
  # order based the sample_id file (NOTE: covDat already ordered based on that)
  id <- unname(sapply(sampleAnn$genoSample_ID, function(x) which(x == colnames(genDat))))
  genDat <- genDat[, id]
  
  # info genotype to match for subset case
  genInfo <- read.table(sprintf('%s%s.txt', genoInfo_file, chr), 
                       header  = T, check.names=F)
  genInfo$match_id <- paste(genInfo$POS, genInfo$REF, genInfo$ALT, sep = '_')
  
  genInfo_model <- read.table(sprintf('%s%s.txt', genoInfo_model_file, chr), 
                        header  = T, check.names=F)
  genInfo_model$match_id <- paste(genInfo_model$POS, genInfo_model$REF, genInfo_model$ALT, sep = '_')
  id_var <- match(genInfo$match_id, genInfo_model$match_id)
  
  geneSnpDist <- readMM(sprintf('%sENSEMBL_gene_SNP_%s_%s_matrix.mtx', InfoFold, cis_ann, chr))
  geneSnpDist <- geneSnpDist[id_var, ]
  ind_SNPs <- lapply(1:ncol(geneSnpDist), function(X) which(geneSnpDist[,X]!=0))
  beta_chr <- betaCoeff[[i]]
  beta_chr <- beta_chr[id_var, ]
  
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
if(!no_zip){system(paste("gzip",sprintf('%spredictedExpression.txt', outFold)))}


