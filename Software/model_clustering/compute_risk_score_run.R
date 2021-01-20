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
parser$add_argument("--type_data", type = "character", help = "")
parser$add_argument("--cases_only", type = "logical", default = FALSE, help = "")
parser$add_argument("--scale_rs", type = "logical", default = FALSE, help = "")
parser$add_argument("--pheno_class_name", type = "character", nargs = '*', help = "")
parser$add_argument("--pvalresFile", type = "character", nargs = '*', help = "file with pvalue results")
parser$add_argument("--corr_thr", type = "double", default = 1, help = "correlation among features threshold")
parser$add_argument("--corrFile", type = "character",  help = "")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
inputFile <- args$inputFile
sampleAnn_file <- args$sampleAnn_file
functR <- args$functR
split_tot <- args$split_tot
pvalresFile <- args$pvalresFile
corr_thr <- args$corr_thr
corrFile <- args$corrFile
type_data <- args$type_data
cases_only <- args$cases_only
min_genes_path <- args$min_genes_path
pheno_class_name <- args$pheno_class_name
scale_rs <- args$scale_rs
outFold <- args$outFold

###################################################################################################################
# sampleAnn_file <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_HARD_clustering/covariateMatrix_CADHARD_All.txt'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_clustering/clustering_functions.R'
# outFold <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Heart_Left_Ventricle/200kb/CAD_GWAS_bin5e-2/German5/devgeno0.01_testdevgeno0/'
# split_tot = 100
# inputFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Heart_Left_Ventricle/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/predictedTscores_splitGenes'
# corrFile <- '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Heart_Left_Ventricle/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/correlation_estimate_tscore.RData'
# pvalresFile <- c('/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Heart_Left_Ventricle/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_Blood_biochemistry_withMed_pheno_covCorr.RData',
# '/psycl/g/mpsziller/lucia/CAD_UKBB/eQTL_PROJECT/OUTPUT_GTEx/predict_CAD/Heart_Left_Ventricle/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_Blood_count_withMed_pheno_covCorr.RData')
# type_data <- 'tscore'
# cases_only <- T
# min_genes_path <- 2
# pheno_class_name <- c('Blood_biochemistry', 'Blood_count')
# corr_thr <- 0.5
# scale_rs <- T
###################################################################################################################

source(functR)

sampleAnn <- read.table(sampleAnn_file, h=T, stringsAsFactors = F)
if(cases_only){
  sampleAnn <- sampleAnn[sampleAnn$Dx == 1,]
}

# note: do not scale internally (treat each cohort seperately)
# load pval res (needed to match features in correlation)
res_pval <- get(load(pvalresFile[1]))
if(type_data == 'tscore'){
  res_pval <- res_pval$tscore
  id_pval <- 8
  id_info <- 2
  id_geno_summ <- 3
}else{
  if(type_data == 'path_Reactome'){
    res_pval <- res_pval$pathScore_reactome
    id_pval <- 13
    id_info <- 1
    id_geno_summ <- 4
  }else{
    if(type_data == 'path_GO'){
      res_pval <- res_pval$pathScore_GO
      id_pval <- 15
      id_info <- 1
      id_geno_summ <- 6
    }else{
      stop('unknown pathway called')
    }
  }
}

# # recompute pvalue if ngenes_tscore > 1
# if(min_genes_path > 1 & grepl('path',type_data)){
#   res_pval <- lapply(res_pval, function(x) x[x$ngenes_tscore >= min_genes_path, ])
#   for(i in 1:length(res_pval)){
#     res_pval[[i]][,id_pval+1] <- qvalue(res_pval[[i]][,id_pval])$qvalues
#     res_pval[[i]][,id_pval+2] <- p.adjust(res_pval[[i]][,id_pval], method = 'BH')
#   }
# }

# load input matrix (new data, old data input saved in clustering result)
if(split_tot == 0){
  
  if(grepl('.txt', inputFile, fixed = TRUE)){
    
    scoreMat <- read.delim(inputFile, h=T, stringsAsFactors = F, check.names = F)
    id_dup <- names(which(table(scoreMat[,1]) > 1)) 
    scoreMat <- scoreMat[!scoreMat[,1] %in% id_dup, ]
    # common_feat <- intersect(scoreMat[,1], res_pval[[1]][, id_info])
    # res_pval <- lapply(res_pval, function(x) x[match(common_feat, x[,id_info]),])
    # scoreMat <- scoreMat[match(common_feat, scoreMat[,1]), ]
    
    rownames(scoreMat) <- scoreMat[,1]
    scoreMat <- scoreMat[,-1]
    if(type_data == 'tscore'){
      new_id <- unname(sapply(colnames(scoreMat), function(x) strsplit(x, split = ' vs reference')[[1]][1]))
      colnames(scoreMat) <- new_id  
    }
    
  }else{
    
    scoreMat <- get(load(inputFile))
    # common_feat <- intersect(scoreMat[,1], res_pval[[1]][, id_info])
    # res_pval <- lapply(res_pval, function(x) x[match(common_feat, x[,id_info]),])
    # scoreMat <- scoreMat[match(common_feat,scoreMat[,1]), ]
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
  # common_feat <- intersect(colnames(scoreMat), res_pval[[1]][, id_info])
  # res_pval <- lapply(res_pval, function(x) x[match(common_feat, x[,id_info]),])
  # scoreMat <- scoreMat[,match(common_feat, colnames(scoreMat))]
  rownames(scoreMat) <- samplesID[[1]]
  # remove sample that have NAs
  id_s <- rowSums(is.na(scoreMat)) == 0
  if(!all(id_s)){scoreMat <- scoreMat[id_s, ]}
  
  common_samples <- intersect(sampleAnn$Individual_ID, rownames(scoreMat))
  sampleAnn <- sampleAnn[match(common_samples, sampleAnn$Individual_ID),]
  scoreMat <- scoreMat[match(common_samples,rownames(scoreMat)),]
  
}
# print(identical(colnames(scoreMat), res_pval[[1]][, id_info]))
print(identical(sampleAnn$Individual_ID, rownames(scoreMat)))

#### filter features based on estimated correlation #####
if(corr_thr < 1){
  
  corr_feat <- get(load(corrFile)) # corrFile estimated based on res_pval, otherwise more intersection/change needed
  corr_feat <- corr_feat$cor
  if(!all(res_pval[[1]][,id_info] %in% rownames(corr_feat))){
    print('not all features in association also in correlation features matrix, possible duplicates')
    res_pval[[1]] <- res_pval[[1]][res_pval[[1]][,id_info] %in% rownames(corr_feat), ]
  }
  common_feat <- intersect(colnames(scoreMat), rownames(corr_feat))
  # clumping: sort according dev_geno (from association results)
  feat_info <- res_pval[[1]][match(common_feat,res_pval[[1]][, id_info]),c(id_info, id_geno_summ)]
  feat_info <- feat_info[order(feat_info[,2], decreasing = T), ]
  corr_feat <- corr_feat[match(feat_info[,1], rownames(corr_feat)), match(feat_info[,1], colnames(corr_feat))]
  
  element_rm <- c()
  stop_cond <- F
  list_feat <- feat_info[,1]
  
  while(!stop_cond){
    
    id <- which(abs(corr_feat[,list_feat[1]])>corr_thr)
    if(length(id)>1){
      element_rm <- c(element_rm, names(id)[-which.max(feat_info[match(names(id),feat_info[,1]), 2])])
    }
    list_feat <- list_feat[!list_feat %in% names(id)]
    # print(length(list_feat))
    stop_cond <- length(list_feat) == 0
    
  }
  
  # for(i in 1:(nrow(corr_feat)-1)){
  #   id <-  which(abs(corr_feat[i:nrow(corr_feat),i])>corr_thr)
  #   if(length(id)>1){
  #     element_rm <- c(element_rm, names(id)[-which.max(res_pval[[1]]$dev_geno[match(names(id), res_pval[[1]][,id_info])])])
  #   }
  # }
  element_rm <- unique(element_rm)
  scoreMat <- scoreMat[, match(common_feat, colnames(scoreMat))] 
  input_data <- scoreMat[, !colnames(scoreMat) %in% element_rm]
  
}else{
  
  input_data <- scoreMat
  
}

###################################################################################################################
# compute PRS
sampleAnn_score <- data.frame(Individual_ID = rownames(input_data))
risk_score <- list()

for(j in 1:length(pvalresFile)){
  
  print(paste( '#########', pheno_class_name[j], '#########'))
  
  res_pval <- get(load(pvalresFile[j]))
  pheno_info <- res_pval$pheno
  
  if(type_data == 'tscore'){
    res_pval <- res_pval$tscore
  }
  if(type_data == 'path_Reactome'){
      res_pval <- res_pval$pathScore_reactome
  }  
  if(type_data == 'path_GO'){
        res_pval <- res_pval$pathScore_GO
  }
  
  # recompute pvalue if ngenes_tscore > 1
  if(min_genes_path > 1 & grepl('path',type_data)){
    res_pval <- lapply(res_pval, function(x) x[x$ngenes_tscore >= min_genes_path, ])
    for(i in 1:length(res_pval)){
      res_pval[[i]][,id_pval+1] <- qvalue(res_pval[[i]][,id_pval])$qvalues
      res_pval[[i]][,id_pval+2] <- p.adjust(res_pval[[i]][,id_pval], method = 'BH')
    }
  }
  
  if(!all(colnames(input_data) %in% res_pval[[1]][,id_info])){stop('not all features in input data also in association analysis')}
  
  risk_score[[j]] <- matrix(ncol = nrow(pheno_info), nrow = nrow(input_data))
  for(i in 1:length(res_pval)){
    print(pheno_info$pheno_id[i])
    tmp <- sapply(1:ncol(input_data), function(x) input_data[, x]*res_pval[[i]][res_pval[[i]][,id_info] == colnames(input_data)[x], id_pval-1])
    colnames(tmp) <- colnames(input_data)
    risk_score[[j]][,i] <- rowSums(tmp)  
  }
  colnames(risk_score[[j]]) <- pheno_info$pheno_id
  if(scale_rs){
    risk_score[[j]] <- scale(risk_score[[j]])
    attr(risk_score[[j]], "scaled:center") <- NULL
    attr(risk_score[[j]], "scaled:scale") <- NULL
  }

}
risk_score <- do.call(cbind, risk_score)
sampleAnn_score <- cbind(sampleAnn_score, risk_score)

write.table(sampleAnn_score, file = sprintf('%s%s_corrThr%s_risk_score_relatedPhenotypes.txt', outFold, type_data, as.character(corr_thr)), col.names = T, row.names = F, sep = '\t', quote = F)

  
tmp <- res_pval[[1]][match(colnames(input_data), res_pval[[1]][,id_info]),c(id_info, id_geno_summ)]
write.table(tmp, file = sprintf('%s%s_features_risk_score_corrThr%s.txt', outFold, type_data, as.character(corr_thr)), col.names = T, row.names = F, sep = '\t', quote = F)




 
