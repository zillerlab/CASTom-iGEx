# ML: evaluate samples FP-FN based on additional endophenotypes

options(stringsAsFactors=F)
options(max.print=1000)
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(coin))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(e1071))
suppressPackageStartupMessages(library(pROC))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(lmtest))

options(bitmapType = 'cairo', device = 'png')

parser <- ArgumentParser(description="ML using zscaled data, evalate missclassified samples")
parser$add_argument("--MLresFile", type = "character", help = "file to be loaded (ML results)")
parser$add_argument("--phenoInfoFile_other", type = "character", help = "")
parser$add_argument("--phenoFile_other", nargs = '*', type = "character", help = "")
parser$add_argument("--name_other", nargs = '*', type = "character", help = "")
parser$add_argument("--phenoFile", type = "character", help = "file to be loaded, additional endophenotypes")
parser$add_argument("--phenoInfoFile", type = "character", help = "file to be loaded, additional endophenotypes info")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--pval_thr", type = "double", default = 1, help = "threshold to filter features")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
MLresFile <- args$MLresFile
phenoFile <- args$phenoFile
phenoInfoFile <- args$phenoInfoFile
phenoFile_other <- args$phenoFile_other
phenoInfoFile_other <- args$phenoInfoFile_other
sampleAnnFile <- args$sampleAnnFile
type_data <- args$type_data
type_input <- args$type_input
name_other <- args$name_other
pval_thr <- args$pval_thr
outFold <- args$outFold


####################################################################################################################
# MLresFile <- 'OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_subset/tscore_original_SVM_FDRpval0.2.RData'
# phenoFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_clustering/phenotypeMatrix_CADsubset.txt'
# phenoFile_other <- c('INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeMatrix_Diet.txt',
#                      'INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeMatrix_Alcohol_use.txt',
#                      'INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeMatrix_Physical_activity.txt',
#                      'INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeMatrix_Arterial_stiffness.txt',
#                      'INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeMatrix_Education.txt')
# name_other <- c('Diet', 'Alcohol_use', 'Physical_activity', 'Arterial_stiffness', 'Education')
# phenoInfoFile_other <- c('INPUT_DATA_GTEx/CAD/Covariates/UKBB/phenotypeDescription_PHESANTproc_CADrelatedpheno.txt')
# phenoInfoFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_clustering/phenotypeDescription_CADsubset.txt'
# sampleAnnFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_clustering/covariateMatrix_CADsubset.txt'
# pval_thr <- 0.2
# type_data <- 'tscore'
# type_input <- 'original'
####################################################################################################################


#######################################################
compute_reg_endopheno <- function(fmla, type_pheno, mat, name_var){
  
  if(type_pheno == 'CONTINUOUS'){
    res <- glm(fmla, data = mat, family = 'gaussian')
    output <- coef(summary(res))[rownames(coef(summary(res))) == name_var,1:4]
  }
  
  if((type_pheno %in% c('CAT_SINGLE_UNORDERED', 'CAT_SINGLE_BINARY', 'CAT_MUL_BINARY_VAR')) | (type_pheno == 'CAT_ORD' & length(unique(na.omit(mat[, 'pheno']))) == 2)){
    
    if(!all(unique(na.omit(mat[, 'pheno']) %in% c(0,1)))){
      
      min_id <- which(mat[, 'pheno'] == min(mat[,'pheno'], na.rm = T))
      mat[min_id,  'pheno'] <- 0
      max_id <- which(mat[, 'pheno'] == max(mat[,  'pheno'], na.rm = T))
      mat[max_id,  'pheno'] <- 1
      
    }
    
    res <- tryCatch(glm(fmla, data = mat, family = 'binomial'),warning=function(...) NA)
    # if(res$converged & is.null(warnings())){
    #   output <- coef(summary(res))[rownames(coef(summary(res))) == name_var,1:4]
    # }else{
    #   output <- rep(NA, 4)
    # }
    if(is.list(res)){
      output <- coef(summary(res))[rownames(coef(summary(res))) == name_var,1:4]
    }else{
      output <- rep(NA, 4)
    }
  }
  
  if(type_pheno == 'CAT_ORD' & length(unique(na.omit(mat[, 'pheno']))) > 2){
    
    mat$pheno <- factor(mat$pheno)
    output <- rep(NA, 4)
    res <- list()
    res$Hess <- NA
    
    res <- tryCatch(polr(fmla, data = mat, Hess=TRUE),warning=function(...) NA)
    if(is.list(res)){
      if(!any(is.na(res$Hess))){
        ct <- coeftest(res)  
        output <- ct[rownames(ct) == name_var,1:4]
      }
    }
    
  }
  
  if(! type_pheno %in% c('CAT_ORD', 'CAT_SINGLE_UNORDERED', 'CAT_SINGLE_BINARY', 'CAT_MUL_BINARY_VAR', 'CONTINUOUS')){
    output <- 'wrong pheno type annotation'
  }
  
  return(output)
  
}
#######################################################

ML_res <- get(load(MLresFile))

sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)
phenoDat <- fread(phenoFile, h=T, stringsAsFactors = F, data.table = F)
phenoInfo <- fread(phenoInfoFile, h=T, stringsAsFactors = F, data.table = F)
phenoInfo_other <- fread(phenoInfoFile_other, h=T, stringsAsFactors = F, data.table = F)

for(i in 1:length(phenoFile_other)){
  
  tmp <-  fread(phenoFile_other[i], h=T, stringsAsFactors = F, data.table = F)
  tmp <- tmp[match(sampleAnn$Individual_ID, tmp$Individual_ID), ]
  print(identical(tmp$Individual_ID,  sampleAnn$Individual_ID))
  tmp <- tmp[, !colnames(tmp) %in% 'Individual_ID']
  
  id_keep <- colSums(is.na(tmp)) != nrow(tmp)
  tmp <- tmp[,id_keep]
  tmp2 <- phenoInfo_other[match(colnames(tmp),phenoInfo_other$pheno_id), ]
  tmp2$pheno_type <- name_other[i]
  phenoInfo <- rbind(phenoInfo, tmp2)
  
  phenoDat <- cbind(phenoDat, tmp)
  
}

df_config <- list()
for(n in 1:5){
  df_config[[n]] <- data.frame(Individual_ID = rownames(ML_res$SVM$nested_cv[[n]]$pred_val), pheno = ML_res$SVM$nested_cv[[n]]$pred_val$y ,
                    pred = ML_res$SVM$nested_cv[[n]]$pred_val$pred, prob = ML_res$SVM$nested_cv[[n]]$pred_val$prob, stringsAsFactors = F)
  df_config[[n]]$fold <- n
}
df_config <- do.call(rbind, df_config)

df_config <- df_config[match(sampleAnn$Individual_ID, df_config$Individual_ID),]
print(identical(df_config$Individual_ID,  sampleAnn$Individual_ID))
print(identical(phenoDat$Individual_ID,  sampleAnn$Individual_ID))

rownames(phenoDat) <- phenoDat$Individual_ID
phenoDat <- phenoDat[, !colnames(phenoDat) %in% 'Individual_ID']

df_config$class <- NA
df_config$class[df_config$pheno == 0 & df_config$pred == 0] <- 'TN'
df_config$class[df_config$pheno == 1 & df_config$pred == 0] <- 'FN'
df_config$class[df_config$pheno == 0 & df_config$pred == 1] <- 'FP'
df_config$class[df_config$pheno == 1 & df_config$pred == 1] <- 'TP'

print(identical(colnames(phenoDat),phenoInfo$pheno_id))

#### check difference in the classification ####
comp_type <- c('FN_vs_TN', 'FP_vs_TN', 'FN_vs_TP', 'FP_vs_TP', 'FP_vs_FN', 'TP_vs_TN')
test_df <- vector(mode = 'list', length = length(comp_type))
for(i in 1:length(comp_type)){
  
  print(paste('#########', comp_type[i], '#########'))
  # test_df[[i]] <- data.frame(pheno_id = colnames(phenoDat), Field = phenoInfo$Field, Coding_meaning = phenoInfo$Coding_meaning, type = phenoInfo$pheno_type)
  # test_df[[i]]$comp <- comp_type[i]
  # test_df[[i]]$pval <- NA
  # test_df[[i]]$pval_BHcorr <- NA
  # test_df[[i]]$statistic <- NA
  # test_df[[i]]$test_type <- NA
  
  comp_single <- strsplit(comp_type[i], split = '_vs_')[[1]]
  id <- which(df_config$class %in% comp_single)
  
  res_glm <- matrix(nrow = ncol(phenoDat),ncol=4)
  for(j in 1:ncol(phenoDat)){
    
    print(j)
    tmp <- cbind(data.frame(pheno = phenoDat[id,j], class = df_config$class[id]), sampleAnn[id, colnames(sampleAnn) %in% c(paste0('PC', 1:10),'Age', 'Gender')])
    tmp$class <- factor(tmp$class, levels = rev(comp_single))
    fmla  <- as.formula(paste('pheno~class+', paste0(colnames(tmp)[-c(1:2)], collapse = '+')))
    res_glm[j,] <- compute_reg_endopheno(mat = tmp, fmla = fmla, type_pheno = phenoInfo$transformed_type[j], name_var = paste0('class', comp_single[1]))
    
    # if(is.integer(phenoDat[,j])){
    #   
    #   if(length(unique(phenoDat[id,j])) > 1){
    #     tmp <- chisq.test(table(phenoDat[id,j], df_config$class[id]))
    #     test_df[[i]]$pval[j] <- tmp$p.value
    #     test_df[[i]]$statistic[j] <-tmp$statistic
    #     test_df[[i]]$test_type[j] <- 'chisq'
    #   }else{
    #     test_df[[i]]$test_type[j] <- 'none'
    #   }
    # }else{
    #   tmp <- wilcox.test(phenoDat[df_config$class == comp_single[1],j] , phenoDat[df_config$class == comp_single[2],j])
    #   test_df[[i]]$pval[j] <- tmp$p.value
    #   test_df[[i]]$statistic[j] <-  tmp$statistic
    #   test_df[[i]]$test_type[j] <- 'wilcox'
    # }
  }
  
  colnames(res_glm) <- c('beta', 'se_beta', 'z', 'pvalue')
  res_glm <- as.data.frame(res_glm)
  res_glm$type_pheno <- phenoInfo$transformed_type
  
  test_df[[i]] <- cbind(data.frame(pheno_id = phenoInfo$pheno_id, Field = phenoInfo$Field, meaning = phenoInfo$Coding_meaning, type = phenoInfo$pheno_type), res_glm)
  test_df[[i]]$pval_corr <- p.adjust(test_df[[i]]$pvalue, method = 'BH')
  test_df[[i]]$comp <- comp_type[i]
}

tot_test <- do.call('rbind', test_df)

# save results
output <- list(test = tot_test, phenoDat = phenoDat, phenoInfo = phenoInfo, sample = df_config)
save(output, file = sprintf('%s%s_%s_SVM_FDRpval%s_misClass.RData', outFold, type_data, type_input, as.character(pval_thr)))

write.table(output$test, file = sprintf('%s%s_%s_SVM_FDRpval%s_misClass.txt', outFold, type_data, type_input, as.character(pval_thr)), col.names = T, row.names = F, sep = '\t', quote = F)




