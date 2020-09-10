# ML: predict pheno status

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

options(bitmapType = 'cairo', device = 'png')


parser <- ArgumentParser(description="ML using zscaled data")
parser$add_argument("--tissues_name", type = "character", nargs = '*', help = "tissues")
parser$add_argument("--inputFile", type = "character", nargs = '*', help = "file to be loaded (predicted tscore or pathScore)")
parser$add_argument("--resMLFile", type = "character", help = "result ML using pval thr")
parser$add_argument("--sampleAnnFile", type = "character", help = "file with samples to be used")
parser$add_argument("--split_tot", type = "integer", default = 0, help = "if 0 then inpuntFile load alone, otherwise splitted version")
parser$add_argument("--pvalresFile", type = "character", nargs = '*', help = "file with pvalue results")
parser$add_argument("--pval_id", type = "integer", default = 0, help = "id to be used on pvalue file")
parser$add_argument("--corr_thr", type = "double", default = -1, help = "correlation among features threshold")
parser$add_argument("--corr_thr_rm", type = "double", default = -1, help = "correlation among features threshold")
parser$add_argument("--pval_thr", type = "double", default = 1, help = "threshold to filter features")
parser$add_argument("--functR", type = "character", help = "functions to be used")
parser$add_argument("--type_data", type = "character", help = "tscore, path_Reactome or path_GO")
parser$add_argument("--type_input", type = "character", default = 'original', help = "original or zscaled")
parser$add_argument("--min_genes_path", type = "integer", default = 1, help = "minimum number of genes for a pathway, if > 1 recompute corrected pvalues")
parser$add_argument("--ncores", type = "integer", default = 5, help = "cores for parallelization")
parser$add_argument("--outFold", type="character", help = "Output file [basename only]")

args <- parser$parse_args()
tissues_name <- args$tissues_name
resMLFile <- args$resMLFile
pvalresFile <- args$pvalresFile
pval_id <- args$pval_id
inputFile <- args$inputFile
sampleAnnFile <- args$sampleAnnFile
split_tot <- args$split_tot
functR <- args$functR
type_data <- args$type_data
type_input <- args$type_input
corr_thr <- args$corr_thr
pval_thr <- args$pval_thr
corr_thr_rm <- args$corr_thr_rm
min_genes_path <- args$min_genes_path
ncores <- args$ncores
outFold <- args$outFold


# ####################################################################################################################
# tissues_name=c('Adipose_Subcutaneous' ,'Adipose_Visceral_Omentum' ,'Adrenal_Gland' ,'Artery_Aorta', 'Artery_Coronary' ,'Colon_Sigmoid' ,'Colon_Transverse' ,
#          'Heart_Atrial_Appendage' ,'Heart_Left_Ventricle' ,'Liver')
# inputFile <- paste0('OUTPUT_GTEx/predict_CAD/',tissues_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/CAD_subset/predictedTscores_splitGenes')
# sampleAnnFile <- 'INPUT_DATA_GTEx/CAD/Covariates/UKBB/CAD_clustering/covariateMatrix_CADsubset.txt'
# split_tot <- 100
# pvalresFile <- paste0('OUTPUT_GTEx/predict_CAD/',tissues_name,'/200kb/CAD_GWAS_bin5e-2/UKBB/devgeno0.01_testdevgeno0/pval_CAD_pheno_covCorr.RData')
# pval_id <- 1
# min_genes_path <- 2
# type_data <- 'tscore'
# outFold <- 'OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_subset/'
# functR <- '/psycl/g/mpsziller/lucia/priler_project/Software/model_ML/ML_functions.R'
# corr_thr <- 0.9
# corr_thr_rm <- 0.5
# pval_thr <- 0.01
# ncores = 5
# type_input <- 'original'
# resMLFile <- 'OUTPUT_GTEx/predict_CAD/AllTissues/200kb/CAD_GWAS_bin5e-2/UKBB/CAD_subset/tscore_original_SVM_FDRpval0.01.RData'
# ####################################################################################################################


source(functR)

res_ML <- get(load(resMLFile))
sampleAnn <- read.table(sampleAnnFile, h=T, stringsAsFactors = F, check.names = F)

res_pval <- list()
scoreMat <- list()
rm_features <- list()
id_rm <- c()
for(i in 1:length(tissues_name)){
  
  print(tissues_name[i])
  
  # load pval res
  res_pval[[i]] <- get(load(pvalresFile[[i]]))
  if(type_data == 'tscore'){
    res_pval[[i]] <- res_pval[[i]]$tscore[[pval_id]]
    id_pval <- 8
    id_info <- 2
  }else{
    if(type_data == 'path_Reactome'){
      res_pval[[i]] <- res_pval[[i]]$pathScore_reactome[[pval_id]]
      id_pval <- 13
      id_info <- 1
    }else{
      if(type_data == 'path_GO'){
        res_pval[[i]] <- res_pval[[i]]$pathScore_GO[[pval_id]]
        id_pval <- 15
        id_info <- 1
      }else{
        stop('unknown pathway called')
      }
    }
  }
  
  # recompute pvalue if ngenes_tscore > 1
  if(min_genes_path > 1 & grepl('path',type_data)){
    res_pval[[i]] <- res_pval[[i]][res_pval[[i]]$ngenes_tscore >= min_genes_path, ]
    res_pval[[i]][,id_pval+1] <- qvalue(res_pval[[i]][,id_pval])$qvalues
    res_pval[[i]][,id_pval+2] <- p.adjust(res_pval[[i]][,id_pval], method = 'BH')
  }
  
  # load input matrix 
  if(split_tot == 0){
    
    scoreMat[[i]] <- get(load(inputFile[[i]]))
    # filter out based on samples and ids
    scoreMat[[i]] <- scoreMat[[i]][match(res_pval[[i]][, id_info],scoreMat[[i]][,1]), ]
    common_samples <- intersect(sampleAnn$Individual_ID, colnames(scoreMat[[i]]))
    
    scoreMat[[i]] <- t(scoreMat[[i]][,match(common_samples,colnames(scoreMat[[i]]))])
    rownames(scoreMat[[i]]) <- common_samples
    colnames(scoreMat[[i]]) <- res_pval[[i]][, id_info]
    
  }else{
    
    ###### load score Mat #######
    scoreMat_list <- vector(mode = 'list', length = split_tot)
    samplesID <- vector(mode = 'list', length = split_tot)
    elementID <- NULL
    
    for(j in 1:split_tot){
      
      print(j)
      if(file.exists(sprintf('%s%i.RData', inputFile[[i]], j))){
        tmp <- get(load(sprintf('%s%i.RData', inputFile[[i]], j)))
        elementID <- c(elementID,tmp[,1])
        samplesID[[j]] <- intersect(sampleAnn$Individual_ID, colnames(tmp))
        scoreMat_list[[j]] <- t(tmp[,match(samplesID[[j]],colnames(tmp))])
      }else{
        print(sprintf('split %i for tissue %s does not exist', j, tissues_name[i]))
      }  
    }
    
    # check samplesID always the same
    if(!all(table(unlist(samplesID)) == split_tot)){print('ERROR: wrong name annotations')}
    
    scoreMat[[i]] <- do.call(cbind, scoreMat_list)
    colnames(scoreMat[[i]]) <- elementID
    rm(scoreMat_list)
    
    # filter out elements that are repeated twice:
    id_dup <- names(which(table(colnames(scoreMat[[i]])) > 1)) 
    scoreMat[[i]] <- scoreMat[[i]][, !colnames(scoreMat[[i]]) %in% id_dup]
    
    id_el <- intersect(colnames(scoreMat[[i]]), res_pval[[i]][,id_info])
    scoreMat[[i]] <- scoreMat[[i]][, match(id_el, colnames(scoreMat[[i]]))]
    res_pval[[i]] <- res_pval[[i]][match(id_el, res_pval[[i]][,id_info]), ]
    
    rownames(scoreMat[[i]]) <- samplesID[[1]]
    # remove sample that have NAs
    id_s <- rowSums(is.na(scoreMat[[i]])) == 0
    if(!all(id_s)){scoreMat[[i]] <- scoreMat[[i]][id_s, ]}
    
    common_samples <- intersect(sampleAnn$Individual_ID, rownames(scoreMat[[i]]))
    scoreMat[[i]] <- scoreMat[[i]][match(common_samples,rownames(scoreMat[[i]])),]
    
  }
  
  print(identical(colnames(scoreMat[[i]]), res_pval[[i]][, id_info]))
  
  # correlation all features
  cor_score <- cor(scoreMat[[i]])
  # remove higly correlated features
  element_rm <- c()
  for(j in 1:(nrow(cor_score)-1)){
    id <-  which(abs(cor_score[j:nrow(cor_score),j])>corr_thr)
    if(length(id)>1){
      element_rm <- c(element_rm, names(id)[-which.min(res_pval[[i]][match(names(id), res_pval[[i]][,id_info]),id_pval])])
    }
  }
  element_rm <- unique(element_rm)
  print(paste(length(element_rm),'features removed due to high correlation'))
  
  # remove features used in the model and correlated ones (0.5)
  model_feat <- res_ML$feat[res_ML$feat$tissue == tissues_name[i], ]
  if(nrow(model_feat)>0){
    tmp_cor <- cor_score[model_feat[, id_info], , drop=F]
    element_rm_f <- model_feat[, id_info]
    for(j in 1:ncol(tmp_cor)){
      if(any(abs(tmp_cor[,j])>corr_thr_rm)){
        element_rm_f <- c(element_rm_f, colnames(tmp_cor)[j])
      }
    }
    element_rm_f <- unique(element_rm_f)
    print(paste(length(element_rm_f),'features removed due to high correlation with pval thr features'))
  }else{
    element_rm_f <- c()  
  }
  rm_features[[i]] <- unique(c(element_rm, element_rm_f))
}

# filter out features from other tissues that are correlated with pval thr
element_rm <- c()
for(i in 1:length(tissues_name)){
  print(i)
  model_feat <- res_ML$feat[res_ML$feat$tissue == tissues_name[i], ]
  tmp_score <- scoreMat[[i]][,model_feat[,id_info], drop = F]
  if(nrow(model_feat)>0){
    for(j in 1:length(tissues_name)){
      if(i!=j){
        print(j)
        tmp_cor <- apply(tmp_score, 2, function(x) apply(scoreMat[[j]], 2, function(y) cor(y,x)))
        for(l in 1:nrow(tmp_cor)){
          if(any(abs(tmp_cor[l,])>corr_thr_rm)){
            element_rm <- c(element_rm, paste0(rownames(tmp_cor)[l], '_tissue_', tissues_name[j]))
          }
        }
      }
    }
  }
}
element_rm <- unique(element_rm)
print(paste(length(element_rm), 'features removed due to high correlation with pval thr features (across tissues)'))

# remove features from scoreMat and filter (same size)
id_rm <- c()
n_rep <- 10
res_pval_rep <- vector(mode = 'list', length = n_rep)
data_rep <- vector(mode = 'list', length = n_rep)

for(j in 1:n_rep){
  
  res_pval_rep[[j]] <- vector(mode = 'list', length = length(tissues_name))
  data_rep[[j]] <- vector(mode = 'list', length = length(tissues_name))
  
  for(i in 1:length(tissues_name)){
  
    model_feat <- res_ML$feat[res_ML$feat$tissue == tissues_name[i], ]
  
    if(nrow(model_feat)>0){
    
      tmp_feat <- sapply(element_rm[grepl(tissues_name[i], element_rm)], function(x) strsplit(x, split = '_tissue_')[[1]][1])
      rm_feat <- unique(c(tmp_feat, rm_features[[i]]))
      # print(length(rm_feat))
      res_pval_rep[[j]][[i]] <- res_pval[[i]][!res_pval[[i]][, id_info] %in% rm_feat, ]
      # exclude additional features with pvalues <=0.1
      res_pval_rep[[j]][[i]] <- res_pval_rep[[j]][[i]][res_pval_rep[[j]][[i]][, id_pval]>=0.1, ]
      
      set.seed(50+i+j*100)
      res_pval_rep[[j]][[i]] <- res_pval_rep[[j]][[i]][sample(1:nrow(res_pval_rep[[j]][[i]]), nrow(model_feat), replace = F), ]
      data_rep[[j]][[i]] <- scoreMat[[i]][, match(res_pval_rep[[j]][[i]][, id_info], colnames(scoreMat[[i]])), drop = F]
    
    }else{
      id_rm <- c(id_rm, i)
    }
  }
}

if(length(id_rm)>0){
  tissues_name <- tissues_name[-id_rm]
  res_pval_rep <- lapply(res_pval_rep, function(x) x[-id_rm])
  data_rep <- lapply(data_rep, function(x) x[-id_rm])
}

res_pval_tot <- vector(mode = 'list', length = n_rep)
for(j in 1:n_rep){
  res_pval_tot[[j]] <- cbind(do.call(rbind, res_pval_rep[[j]]), data.frame(tissue = unlist(lapply(1:length(tissues_name), function(x) rep(tissues_name[x], nrow(res_pval_rep[[j]][[x]]))))))
  res_pval_tot[[j]]$new_id <- paste(res_pval_tot[[j]][, id_info], 'tissue', res_pval_tot[[j]]$tissue, sep = '_')
}

print('load data: completed')


#############################################################################################################
# for each repetition run SVM
output <- vector(mode = 'list', length = n_rep)
for(j in 1:n_rep){
  
  #### normalize features tot ####
  tot_data <- list()
  par_norm <- list()
  for(i in 1:length(tissues_name)){
    # print(tissues_name[i])
    tmp <- normalize_feat(data_rep[[j]][[i]])
    tot_data[[i]] <- tmp$data
    par_norm[[i]] <- tmp$par
    
    if(type_input == 'zscaled'){
      tot_data[[i]] <- sapply(1:ncol(tot_data[[i]]), function(x) tot_data[[i]][, x]*res_pval_rep[[j]][[i]][res_pval_rep[[j]][[i]][,id_info] == colnames(tot_data[[i]])[x], id_pval-1])
    }
    colnames(tot_data[[i]]) <- paste0(colnames(data_rep[[j]][[i]]), '_tissue_', tissues_name[i])
  }
  print('normalization data: completed')
  
  #### create SVM model ####
  input_data <- as.data.frame(do.call(cbind, tot_data))
  SVM_tot <- SVM_model(data = input_data, info_sample = sampleAnn, n_cv = 5, ncores = ncores)
  print(SVM_tot$performance_nestedCV)
  print('SVM data: completed')
  
  ### save ####
  output[[j]] <- list(feat = res_pval_tot[[j]], SVM_nestedCV = SVM_tot$performance_nestedCV, SVM_CV = SVM_tot$performance_CV)
  
  print('save results: completed')
  
}

res_tot <- list(output_rep = output, data = scoreMat)
save(res_tot, file = sprintf('%s%s_%s_SVM_FDRpval%s_randomFeat.RData', outFold, type_data, type_input, as.character(pval_thr)))


