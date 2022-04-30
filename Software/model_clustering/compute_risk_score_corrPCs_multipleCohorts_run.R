# compute gene PRS/path PRS

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
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

parser <- ArgumentParser(description="gene-RS, correct for PCs and across multiple cohorts")
parser$add_argument("--genes_to_filter", type = "character", default = NULL, help = "additional file to filter genes")
parser$add_argument("--sampleAnn_file", type = "character", nargs = '*', help = "")
parser$add_argument("--name_cohorts", type = "character", nargs = '*', help = "")
parser$add_argument("--inputFile", type = "character", nargs = '*', help = "input files (scores)")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "")
parser$add_argument("--cases_only", type = "logical", default = FALSE, help = "")
parser$add_argument("--scale_rs", type = "logical", default = FALSE, help = "")
parser$add_argument("--pheno_class_name", type = "character", nargs = '*', help = "")
parser$add_argument("--pvalresFile", type = "character", nargs = '*', help = "file with pvalue results")
parser$add_argument("--sqcorr_thr", type = "double", default = 1, help = "squared correlation among features threshold")
parser$add_argument("--n_max", type = "integer", default = 50000, help = "maximum number of samples to split data")
parser$add_argument("--corrFile", type = "character",  help = "")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
name_cohorts <- args$name_cohorts
genes_to_filter <- args$genes_to_filter
inputFile <- args$inputFile
sampleAnn_file <- args$sampleAnn_file
functR <- args$functR
split_tot <- args$split_tot
pvalresFile <- args$pvalresFile
sqcorr_thr <- args$sqcorr_thr
corrFile <- args$corrFile
type_data <- args$type_data
cases_only <- args$cases_only
pheno_class_name <- args$pheno_class_name
scale_rs <- args$scale_rs
n_max <- args$n_max
min_genes_path <- args$min_genes_path
outFold <- args$outFold

###################################################################################################################
# genes_to_filter='compare_prediction_UKBB_SCZ-PGC/DLPC_CMC_filter_genes_matched_datasets.txt'
# name_cohorts <- read.table('/home/luciat/eQTL_PROJECT/INPUT_DATA/SCZ_cohort_names_CLUST', header = F, stringsAsFactors = F)$V1[c(1:2, 16)]
# sampleAnn_file=paste0('/home/luciat/eQTL_PROJECT/INPUT_DATA/Covariates/',name_cohorts,'.covariateMatrix_old.txt')
# functR <- '/home/luciat/castom-igex/Software/model_clustering/clustering_functions.R'
# split_tot = 0
# pheno_class_name = read.table('/home/luciat/UKBB_SCZrelated/name_pheno_risk_score.txt', h=F,stringsAsFactors = F)$V1
# inputFile = paste0('/home/luciat/eQTL_PROJECT/OUTPUT_CMC/predict_PGC/200kb/',name_cohorts, '/devgeno0.01_testdevgeno0/predictedTscores.txt')
# corrFile = '/home/luciat/UKBB_SCZrelated/DLPC_CMC/correlation_estimate_tscore.RData'
# type_data='tscore'
# sqcorr_thr = 0.1
# scale_rs <- T
# min_genes_path <- 2
# pvalresFile=sprintf('/home/luciat/UKBB_SCZrelated/DLPC_CMC/pval_%s_pheno_covCorr.RData', pheno_class_name )
# cases_only=T
# n_max = 50000
###################################################################################################################

source(functR)

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

if(!is.null(genes_to_filter)){
  genes_filt <- read.table(genes_to_filter, h=T, stringsAsFactors = F, sep = '\t')
  genes_filt <- genes_filt[genes_filt$keep & !is.na(genes_filt$keep),]
  # the first one is used to filter input data
  res_pval[[1]] <- res_pval[[1]][res_pval[[1]]$ensembl_gene_id %in% genes_filt$ensembl_gene_id,]
}


### load score data ###
sampleAnn <- vector(mode = 'list', length = length(name_cohorts))
scoreMat <- vector(mode = 'list', length = length(name_cohorts))

for(c_id in 1:length(name_cohorts)){
  
  print(name_cohorts[c_id])
  
  sampleAnn[[c_id]] <- read.table(sampleAnn_file[[c_id]], h = T, stringsAsFactors = F)
  if(cases_only){
    sampleAnn[[c_id]] <- sampleAnn[[c_id]][sampleAnn[[c_id]]$Dx == 1,]
  }
  
  sampleAnn[[c_id]]$Temp_ID <- sampleAnn[[c_id]]$Individual_ID
  sampleAnn[[c_id]]$cohort <- rep(name_cohorts[c_id], nrow(sampleAnn[[c_id]]))
  
  load_output <- load_input_matrix(inputFile[c_id], sampleAnn[[c_id]], res_pval[[1]], split_tot, id_info)
  scoreMat[[c_id]] <- load_output$scoreMat
  sampleAnn[[c_id]] <- load_output$sampleAnn
  
  # remove sample that have NAs
  id_s <- rowSums(is.na(scoreMat[[c_id]])) == 0
  if(!all(id_s)){scoreMat[[c_id]] <- scoreMat[[c_id]][id_s,]}
  sampleAnn[[c_id]] <- sampleAnn[[c_id]][match(rownames(scoreMat[[c_id]]), 
                                               sampleAnn[[c_id]]$Individual_ID), ]
  
  print(dim(scoreMat[[c_id]]))
  print(dim(sampleAnn[[c_id]]))
  
}
try(if(any(!sapply(scoreMat, function(x) all(colnames(x) %in% res_pval[[1]][, id_info])))) 
  stop("ERROR: wrong filtering elements"))
res_pval[[1]] <- res_pval[[1]][match(colnames(scoreMat[[1]]), res_pval[[1]][, id_info]),]
print(identical(colnames(scoreMat[[1]]), res_pval[[1]][,id_info]))

scoreMat_tot <- do.call(rbind, scoreMat)
sampleAnn_list <- sampleAnn
sampleAnn <- do.call(rbind, sampleAnn)
# match to have the same samples and same order with annotation
sampleAnn <- sampleAnn[match(rownames(scoreMat_tot), sampleAnn$Temp_ID), ]
sampleAnn$cohort_id <- as.numeric(as.factor(sampleAnn$cohort))

if(sqcorr_thr < 1){
  
  corr_feat <- get(load(corrFile)) # corrFile estimated based on res_pval, otherwise more intersection/change needed
  corr_feat <- corr_feat$cor
  if(!all(res_pval[[1]][,id_info] %in% rownames(corr_feat))){
    print('not all features in association also in correlation features matrix, possible duplicates')
    res_pval[[1]] <- res_pval[[1]][res_pval[[1]][,id_info] %in% rownames(corr_feat), ]
  }
  common_feat <- intersect(colnames(scoreMat_tot), rownames(corr_feat))
  common_feat <- intersect(common_feat, res_pval[[1]][,id_info])
  # clumping: sort according dev_geno (from association results)
  tmp <- res_pval[[1]][match(common_feat,res_pval[[1]][, id_info]),]
  
  element_rm <- clumping_features(res_pval = tmp, id_info = id_info, 
                                  corr_feat = corr_feat, 
                                  id_pval = id_geno_summ, 
                                  corr_thr = sqrt(sqcorr_thr))
  
  scoreMat_tot <- scoreMat_tot[, match(common_feat, colnames(scoreMat_tot))] 
  input_data_notcorr <- scoreMat_tot[, !colnames(scoreMat_tot) %in% element_rm]
  
}else{
  input_data_notcorr <- scoreMat_tot
}
rm(scoreMat_tot)

## correct for PCs ##
# remove PCs1-10 for each genes
input_data <- matrix(ncol = ncol(input_data_notcorr), nrow = nrow(input_data_notcorr))
rownames(input_data) <- rownames(input_data_notcorr)
colnames(input_data) <- colnames(input_data_notcorr)

name_cov <- setdiff(colnames(sampleAnn),
                    c('Individual_ID', 'genoSample_ID', 'Dx', 'Sex',
                      'Age', 'Temp_ID', 'cohort', 'cohort_id'))
fmla <- as.formula(paste('g ~', paste0(name_cov, collapse = '+')))

for(i in 1:ncol(input_data_notcorr)){
  # print(i)
  tmp <- data.frame(g = input_data_notcorr[,i], sampleAnn[, name_cov])
  reg <- lm(fmla, data = tmp)
  input_data[,i] <- reg$residuals
}
print("corrected for PCs")


### compute PRS ##
sampleAnn_score <- data.frame(Individual_ID = rownames(input_data))
risk_score <- list()

for(j in 1:length(pvalresFile)){
  
  print(paste( '#########', pheno_class_name[j], '#########'))
  
  res_pval <- get(load(pvalresFile[j]))
  rm(final)
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
  
  if(!all(colnames(input_data) %in% res_pval[[1]][,id_info])){
    stop('not all features in input data also in association analysis')}
  
  risk_score[[j]] <- matrix(ncol = nrow(pheno_info), nrow = nrow(input_data))
  
  for(i in 1:length(res_pval)){
    print(pheno_info$pheno_id[i])
    
    tmp_z <- res_pval[[i]][match(colnames(input_data), res_pval[[i]][,id_info]),id_pval-1]
    
    if(nrow(input_data) > n_max){
      
      n_split <- floor(nrow(input_data)/100)
      split_row <- lapply(1:100, function(x) (1:n_split) + n_split*(x-1))
      if(max(split_row[[100]]) < nrow(input_data)){
        split_row[[100]] <- c(split_row[[100]], (split_row[[100]][n_split]+1):nrow(input_data))
      }
      if(max(split_row[[100]]) > nrow(input_data)){
        split_row[[100]] <- split_row[[100]][split_row[[100]] %in% 1:nrow(input_data)]
      }
      tmp <- list()
      for(l in 1:length(split_row)){
        tmp[[l]] <- input_data[split_row[[l]],]*matrix(rep(tmp_z, length(split_row[[l]])), byrow = T, nrow = length(split_row[[l]]))
      }
      tmp <- do.call(rbind, tmp)
      colnames(tmp) <- colnames(input_data)
      
    }else{
      # tmp <- sapply(1:ncol(input_data), function(x) input_data[, x]*res_pval[[i]][res_pval[[i]][,id_info] == colnames(input_data)[x], id_pval-1])
      # colnames(tmp) <- colnames(input_data)
      tmp <- input_data*matrix(rep(tmp_z, nrow(input_data)), byrow = T, nrow = nrow(input_data))
      colnames(tmp) <- colnames(input_data)
    }
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

write.table(sampleAnn_score, file = sprintf('%s%s_corr2Thr%s_risk_score_relatedPhenotypes.txt', outFold, type_data, as.character(sqcorr_thr)), col.names = T, row.names = F, sep = '\t', quote = F)

tmp <- res_pval[[1]][match(colnames(input_data), res_pval[[1]][,id_info]),c(id_info, id_geno_summ)]
write.table(tmp, file = sprintf('%s%s_features_risk_score_corr2Thr%s.txt', outFold, type_data, as.character(sqcorr_thr)), col.names = T, row.names = F, sep = '\t', quote = F)




