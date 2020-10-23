# compute gene PRS/path PRS

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SparseM))
suppressPackageStartupMessages(library(SNFtool))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(ggplot2))
options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="gene-RS/path-RS for ax external cohort")
parser$add_argument("--sampleAnn_file", type = "character", help = "")
parser$add_argument("--inputFile", type = "character", help = "input files (scores)")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--pvalresFile", type = "character", default = 'NA', help = "file with pvalue results")
parser$add_argument("--pval_id", type = "integer", default = 0, help = "id to be used on pvalue file")
parser$add_argument("--corr_thr", type = "double", default = -1, help = "correlation among features threshold")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFile <- args$inputFile
sampleAnn_file <- args$sampleAnn_file
functR <- args$functR
split_tot <- args$split_tot
type_cluster <- args$type_cluster
pvalresFile <- args$pvalresFile
corr_thr <- args$corr_thr
pval_id <- args$pval_id
type_data <- args$type_data
min_genes_path <- args$min_genes_path
outFold <- args$outFold

# ###################################################################################################################
# sampleAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/German5/covariateMatrix.txt'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Heart_Left_Ventricle/200kb/CAD_GWAS_bin5e-2/German5/devgeno0.01_testdevgeno0/'
# split_tot = 0
# inputFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Heart_Left_Ventricle/200kb/CAD_GWAS_bin5e-2/German5/devgeno0.01_testdevgeno0/predictedTscores.txt'
# pvalresFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Heart_Left_Ventricle/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData'
# type_data <- 'tscore'
# pval_id <- 1
# corr_thr <- 0.9
# min_genes_path <- 2
# ###################################################################################################################

source(functR)

sampleAnn <- read.table(sampleAnn_file, h=T, stringsAsFactors = F)

# note: do no scale internally (treat each cohort seperately)

# load pval res
res_pval <- get(load(pvalresFile))
if(type_data == 'tscore'){
  res_pval <- res_pval$tscore[[pval_id]]
  id_pval <- 8
  id_info <- 2
}else{
  if(type_data == 'path_Reactome'){
    res_pval <- res_pval$pathScore_reactome[[pval_id]]
    id_pval <- 13
    id_info <- 1
  }else{
    if(type_data == 'path_GO'){
      res_pval <- res_pval$pathScore_GO[[pval_id]]
      id_pval <- 15
      id_info <- 1
    }else{
      stop('unknown pathway called')
    }
  }
}

# recompute pvalue if ngenes_tscore > 1
if(min_genes_path > 1 & grepl('path',type_data)){
  res_pval <- res_pval[res_pval$ngenes_tscore >= min_genes_path, ]
  res_pval[,id_pval+1] <- qvalue(res_pval[,id_pval])$qvalues
  res_pval[,id_pval+2] <- p.adjust(res_pval[,id_pval], method = 'BH')
}

# load input matrix (new data, old data input saved in clustering result)
if(split_tot == 0){
  
  if(grepl('.txt', inputFile, fixed = TRUE)){
    scoreMat <- read.delim(inputFile, h=T, stringsAsFactors = F, check.names = F)
    id_dup <- names(which(table(scoreMat[,1]) > 1)) 
    scoreMat <- scoreMat[!scoreMat[,1] %in% id_dup, ]
    res_pval <- res_pval[!res_pval[,id_info] %in% id_dup,]
    scoreMat <- scoreMat[match(res_pval[, id_info],scoreMat[,1]), ]
    rownames(scoreMat) <- scoreMat[,1]
    scoreMat <- scoreMat[,-1]
    if(type_data == 'tscore'){
      new_id <- unname(sapply(colnames(scoreMat), function(x) strsplit(x, split = ' vs reference')[[1]][1]))
      colnames(scoreMat) <- new_id  
    }
  }else{
    scoreMat <- get(load(inputFile))  
    scoreMat <- scoreMat[match(res_pval[, id_info],scoreMat[,1]), ]
    rownames(scoreMat) <- scoreMat[,1]
    scoreMat <- scoreMat[,-1]
  }
  
  common_samples <- intersect(sampleAnn$Individual_ID, colnames(scoreMat))
  sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
  scoreMat <- t(scoreMat[,match(common_samples,colnames(scoreMat))])
  
}else{
  
  ###### load score Mat #######
  scoreMat_list <- vector(mode = 'list', length = split_tot)
  samplesID <- vector(mode = 'list', length = split_tot)
  elementID <- NULL
  
  for(i in 1:split_tot){
    
    print(i)
    tmp <- get(load(sprintf('%s%i.RData', inputFile, i)))
    elementID <- c(elementID,tmp[,1])
    samplesID[[i]] <- intersect(sampleAnn$Individual_ID, colnames(tmp))
    scoreMat_list[[i]] <- t(tmp[,match(samplesID[[i]],colnames(tmp))])
    
  }
  
  # check samplesID always the same
  if(!all(table(unlist(samplesID)) == split_tot)){print('ERROR: wrong name annotations')}
  
  scoreMat <- do.call(cbind, scoreMat_list)
  colnames(scoreMat) <- elementID
  rm(scoreMat_list)
  
  # filter out elements that are repeated twice:
  id_dup <- names(which(table(colnames(scoreMat)) > 1)) 
  scoreMat <- scoreMat[, !colnames(scoreMat) %in% id_dup]
  
  scoreMat <- scoreMat[, match(id_el, colnames(scoreMat))]
  
  rownames(scoreMat) <- samplesID[[1]]
  # remove sample that have NAs
  id_s <- rowSums(is.na(scoreMat)) == 0
  if(!all(id_s)){scoreMat <- scoreMat[id_s, ]}
  
  common_samples <- intersect(sampleAnn$Individual_ID, rownames(scoreMat))
  sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
  scoreMat <- scoreMat[match(common_samples,rownames(scoreMat)),]
  
}
print(identical(colnames(scoreMat), res_pval[, id_info]))
print(identical(sampleAnn$Individual_ID, rownames(scoreMat)))

input_data <- scoreMat
tmp <- sapply(1:ncol(input_data), function(x) input_data[, x]*res_pval[res_pval[,id_info] == colnames(input_data)[x], id_pval-1])
colnames(tmp) <- colnames(scoreMat)
risk_score <- rowSums(tmp)

sampleAnn_score <- data.frame(Individual_ID = names(risk_score), Dx = sampleAnn$Dx, score = risk_score)

write.table(sampleAnn_score, file = sprintf('%s%s_risk_score_model%s.txt', outFold, type_data, model_name), col.names = T, row.names = F, sep = '\t', quote = F)
